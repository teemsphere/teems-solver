# TEEMS Solver — Image Building and Publishing Guide

How to build and publish the solver Docker images across architectures
and ISA performance levels. Written for the maintainer who publishes
`matthewcantele/teems_base`; the expedited/full builds are what end
users run locally (the HSL licences forbid distributing the final
image, so only the **base** image is ever published).

## 1. The image family

| Image | Built by | Distributed? | Contents |
|---|---|---|---|
| `matthewcantele/teems_base` | maintainer, `docker/base_build/Dockerfile` | **yes** (Docker Hub) | Debian bookworm + toolchain + MPICH + PETSc + a `src/` snapshot |
| `teems:<tag>` (expedited) | end user, `docker/expedited_build/Dockerfile` | no (HSL licence) | base image + user's HSL tarballs + compiled solver, squashed to a slim runtime |
| `teems:<tag>` (full) | end user, `docker/full_build/Dockerfile` | no | same result without the pre-built base (compiles everything, ~40 min) |

teems-R invokes `teems:<docker_tag>` (default `latest`; override with
`ems_option_set(docker_tag = ...)`) and runs
`/opt/teems-solver/solver/hsl` inside it.

## 2. ISA levels: how performance is selected

Every layer is compiled at one **ISA level** controlled by the `MARCH`
build argument. The base image resolves `MARCH` to concrete compiler
flags, uses them for MPICH/PETSc, and records them in
`/opt/teems-solver/archflags`. The expedited build reads that file and
compiles the HSL libraries and the solver at the *same* level — the
two layers cannot mismatch (the pre-2026-07 images had exactly that
bug: base libs at `-march=native`, solver objects at plain baseline).

| `MARCH` value | Arch | Runs on | Notes |
|---|---|---|---|
| `x86-64-v2` (amd64 default) | amd64 | every Intel/AMD since ~2009 (Nehalem/Bulldozer: SSE4.2, popcnt) | the safe published default |
| `x86-64-v3` | amd64 | Intel Haswell 2013+ / AMD Zen 1+ (AVX2, FMA) | the recommended performance tag; most user machines qualify |
| `x86-64-v4` | amd64 | AVX-512 chips (Intel server/12th-gen-, AMD Zen 4+) | publish only if the A/B gate shows a win over v3 |
| `armv8-a` (arm64 default) | arm64 | every 64-bit ARM incl. all Apple M-series, Graviton | NEON is part of the baseline |
| `armv8.2-a`, `armv8.4-a`, `armv8-a+sve` … | arm64 | newer ARM (M1+ is ≥ v8.4; Graviton3+ has SVE) | benchmark before publishing — see §7 |

Two performance facts to keep in mind when deciding how many variants
to publish:

- **The BLAS3 inner kernels don't need any of this.** MA48's dense
  frontal operations (`dgemm`/`dtrsm`) resolve at runtime to Debian's
  OpenBLAS, which dispatches to CPU-specific kernels dynamically. The
  `MARCH` level affects everything *around* them: MA48's pivot/sparse
  code, HSL_MP48, the solver's own C (`formula_eval`,
  `jacobian_fill`), and PETSc assembly.
- **arm64 levels buy less than amd64 levels.** The amd64 ladder widens
  vectors (SSE→AVX2→AVX-512); on arm64, 128-bit NEON is already in the
  baseline and Apple M-series has no SVE, so `armv8.4-a` mostly adds
  minor instructions. Expect the measurable arm64 win to be small;
  publish extra arm64 tags only with A/B evidence (Graviton3+/SVE is
  the case most likely to pay).

## 3. One-time setup

```bash
# Docker Hub account with push rights to matthewcantele/teems_base
docker login

# a buildx builder (docker-container driver; required for --push and
# multi-arch manifests). One already exists on the dev box:
docker buildx ls           # look for teems-builder
# if absent:
docker buildx create --name teems-builder --driver docker-container --bootstrap
```

Cross-building the arm64 base from an amd64 host goes through QEMU
emulation: it works, but MPICH + PETSc compile for **hours** and
emulation has historically been the flaky path. Strongly prefer a
native arm64 builder when one is available — an Apple-silicon Mac
running Docker Desktop, an AWS Graviton instance, or GitHub Actions'
free arm64 runners — and combine the per-arch pushes into one manifest
afterwards (§5). QEMU remains acceptable for occasional rebuilds if no
native machine is at hand.

## 4. Building the base variants

Run from the repository root (the Dockerfile `COPY ./src` needs the
repo as context). One build per (platform, MARCH) pair:

```bash
V=1.1.0   # base image version — bump on every published change

# amd64, portable default (x86-64-v2)
docker buildx build --builder teems-builder \
  --platform linux/amd64 \
  -t matthewcantele/teems_base:x86-64-v2 \
  -t matthewcantele/teems_base:${V}-x86-64-v2 \
  -f docker/base_build/Dockerfile --push .

# amd64, performance tier (AVX2/FMA)
docker buildx build --builder teems-builder \
  --platform linux/amd64 --build-arg MARCH=x86-64-v3 \
  -t matthewcantele/teems_base:x86-64-v3 \
  -t matthewcantele/teems_base:${V}-x86-64-v3 \
  -f docker/base_build/Dockerfile --push .

# arm64, portable default (covers all M-series and Graviton)
docker buildx build --builder teems-builder \
  --platform linux/arm64 \
  -t matthewcantele/teems_base:armv8-a \
  -t matthewcantele/teems_base:${V}-armv8-a \
  -f docker/base_build/Dockerfile --push .
```

Notes:

- `--build-arg MARCH=...` overrides the per-arch default
  (`x86-64-v2` / `armv8-a`). Never pass `native` for a published
  image — that bakes in the build host's CPU.
- Each variant takes ~40 min natively (MPICH + PETSc from source).
- The version-pinned tags (`${V}-…`) are immutable history; the plain
  level tags are what `latest` and users track.
- `-march` changes do **not** preserve bit-identical solver output:
  plan a golden re-anchoring (§6) whenever the published levels
  change.

## 5. Publishing: the `latest` manifest

`latest` should be a multi-arch manifest whose per-arch entries are
the *portable* variants, so a bare `FROM matthewcantele/teems_base`
is always safe:

```bash
docker buildx imagetools create -t matthewcantele/teems_base:latest \
  matthewcantele/teems_base:x86-64-v2 \
  matthewcantele/teems_base:armv8-a

# verify what got published
docker buildx imagetools inspect matthewcantele/teems_base:latest
```

Publishing checklist:

1. All variant builds pushed (§4) and the `latest` manifest updated.
2. `docker pull` each tag on a clean machine (or after
   `docker image rm`) and confirm `docker run --rm <tag> cat
   /opt/teems-solver/archflags` prints the expected level.
3. Build the expedited image from each published base variant and run
   the golden suite against it (§6).
4. On the oldest amd64 hardware available, run one expedited solve
   from the `x86-64-v2` base — this is the SIGILL regression test.
5. Bump the version labels in `docker/base_build/Dockerfile`
   (`org.opencontainers.image.version`) and
   `docker/base_build/README.md`; update the manual's dependencies
   page if user-facing commands changed.

## 6. The expedited build (what users — and the validation flow — run)

```bash
# from the repo root, with the four HSL tarballs in hsl/
docker build --pull \
  --build-arg BASE_IMAGE=matthewcantele/teems_base:x86-64-v3 \
  -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
  -f ./docker/expedited_build/Dockerfile \
  .
```

- `BASE_IMAGE` selects the variant; omit it for the portable
  `:latest` base. The solver and HSL libraries are automatically
  compiled at the level recorded in the base's `archflags` (a base
  image predating that file falls back to level-less flags — old
  behaviour, no error).
- Building against an old (pre-archflags) base is supported but
  reproduces the historical layer mismatch; rebuild the base first
  when performance matters.
- To hold several variants side by side, tag them apart
  (`teems:v3`, `teems:v2`) and point teems-R at one with
  `ems_option_set(docker_tag = "v3")`.

**Validation** (maintainer, after any flag or dependency change):

1. Rebuild the `teems-audit` container from the new base and run
   `.audit/verify.sh`. Flag changes are *expected* to break
   bit-identity — re-anchor the 9 golden manifests once per published
   flag scheme, not per variant build.
2. Run the teems-R test suite against the new `teems:latest`
   (`test-ems_solve.R` exercises deploy → solve → compose end to
   end).
3. A/B the performance claim that justified the variant (§7).

## 7. Benchmark gates

Per the roadmap's standing rule, a performance-motivated variant needs
interleaved A/B evidence on an idle machine (`.audit/bench_run.sh`
rigs; min-of-N, paired runs, discard first-run-of-session walls):

- `x86-64-v3` vs `x86-64-v2` on bench-static and bench-inter-L —
  expected win concentrated in `formula_eval`/assembly, since OpenBLAS
  already dispatches optimally in both.
- `x86-64-v4` vs `v3` before ever publishing a v4 tag (AVX-512
  downclocking can make it a net loss).
- Any arm64 tier vs `armv8-a` (see §2 — burden of proof is on the
  extra tag).
- **OpenBLAS thread pool**: the runtime image ships Debian's
  pthread-variant OpenBLAS with no thread cap. With several MPI ranks
  each opening a BLAS pool this can oversubscribe; A/B
  `OPENBLAS_NUM_THREADS=1` (settable via `docker run -e`) on the
  bench rigs and bake the winner into the expedited runtime stage as
  an `ENV`.

Record results in `docs/solver-reference.md` §9 and flip the
corresponding roadmap 6.6 rows.

## 8. What lives where (quick reference)

| Concern | File | Knob |
|---|---|---|
| base ISA level | `docker/base_build/Dockerfile` | `--build-arg MARCH=` |
| level handoff to expedited | `${BUILD_DIR}/archflags` in the base image | written at base build |
| MPICH production flags | base/full Dockerfiles | `--enable-fast=all,O3`, `MPICHLIB_*FLAGS` |
| PETSc flags | base/full Dockerfiles | `COPTFLAGS`/`FOPTFLAGS` (= `-Ofast` + arch), `--with-memalign=64` |
| solver objects | expedited/full Dockerfiles | `make OPT="-Ofast $ARCH_FLAGS"` |
| HSL libraries | expedited/full Dockerfiles | `configure CFLAGS/FFLAGS="-O2 $ARCH_FLAGS"` (deliberately not `-Ofast`: keep MA48's pivoting semantics conservative) |
| base variant used by expedited | `docker/expedited_build/Dockerfile` | `--build-arg BASE_IMAGE=` |
| image teems-R runs | R option | `ems_option_set(docker_tag = ...)` |

Rebuild + republish the **base** when: Debian base, MPICH/PETSc
versions, or the flag scheme change. Solver `src/` changes do *not*
require a base republish — the expedited build overwrites the baked
`src/` snapshot with the current tree.

Pinned versions: MPICH 5.0.1, PETSc v3.25.0 (build args in the
Dockerfiles); upstream sources for citation live in `teems-docs/mpich`
and `teems-docs/petsc`.
