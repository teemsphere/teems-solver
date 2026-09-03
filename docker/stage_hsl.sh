#!/bin/bash
# stage_hsl.sh <libHSL-tarball> <staging-root>
#
# Stage the HSL sources the solver build consumes from a libHSL source
# snapshot (single tarball from https://licences.stfc.ac.uk), producing
# the lib/<pkg>/... layout the Dockerfiles compile from:
#
#   lib/ma48/src/{ma48d.f,ddeps.f}            -> libma48.so (dynamic)
#   lib/ma51/src/{ma51d.f,ddeps.f}            -> libma51.so (dynamic)
#   lib/mp48/src/{hsl_mp48d.f90,ddeps.f,common90.f90}
#                                             -> staged into solver src
#                                                (TEEMS patches applied
#                                                there by mp48_mod.sh);
#                                                common90 -> hsl_mp01.mod
#   lib/mc66/src/{hsl_mc66d.f90,ddeps90.f90}  -> staged into solver src
#   lib/mc79/src/{hsl_mc79i.f90,hsl_mc79i_ciface.f90}
#   lib/mc79/include/hsl_mc79i.h              -> staged into solver src
#   lib/ma60/src/ma60d.f, lib/mc71/src/mc71d.f, lib/fd15/fd15d.f
#                                             -> -condest diagnostics
#   LICENSES/hsl/LICENCE, LICENSES/hsl/VERSION
#
# The ddeps files bundled with the historical per-package tarballs are
# exact concatenations of the individual dependency decks; the recipes
# below reproduce them byte-for-byte from libHSL 2026.8.4 (except
# mc66's ddeps90, which differs only in module-name case). Every input
# deck is checksum-pinned: a future libHSL snapshot that changes any
# used deck fails loudly here and must be re-verified (diff the decks,
# re-run the parity kits) before its checksums are recorded. The MP48
# sources additionally carry TEEMS patches applied by src/mp48_mod.sh,
# which also fails loudly on drift.
set -euo pipefail

TARBALL=${1:?usage: stage_hsl.sh <libHSL-tarball> <staging-root>}
STAGE=${2:?usage: stage_hsl.sh <libHSL-tarball> <staging-root>}

if [ ! -f "$TARBALL" ]; then
    echo "stage_hsl.sh: '$TARBALL' not found." >&2
    echo "Download the libHSL source tarball (libHSL.v<version>.tar.gz) from" >&2
    echo "https://licences.stfc.ac.uk and place it in the hsl/ directory." >&2
    echo "Individual per-package HSL tarballs are no longer a supported input." >&2
    exit 1
fi
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
tar xzf "$TARBALL" -C "$TMP"
ROOT=$(echo "$TMP"/libHSL.v*)
if [ ! -d "$ROOT" ]; then
    echo "stage_hsl.sh: '$TARBALL' is not a libHSL source tarball (no" >&2
    echo "libHSL.v<version>/ directory inside). Individual per-package HSL" >&2
    echo "tarballs are no longer a supported input: download the libHSL" >&2
    echo "source snapshot from https://licences.stfc.ac.uk instead; see" >&2
    echo "README (HSL Libraries)." >&2
    exit 1
fi
VERSION=$(basename "$ROOT" | sed 's/^libHSL\.v//')

# Checksums of every deck consumed, as verified against the historical
# per-package tarballs (byte-identical) for libHSL 2026.8.4.
MANIFEST="\
96cf828b2b77b3087ac82c4c8c31d706b0f0de1d41da5863183f25bd015c470c  ma48/ma48d.f
ed2a4c5f31422272ad908a71dfe9a99bf2d8e273a128da5fd421c25bcac018e7  ma51/ma51d.f
0c0fe7d928ca96e18d49a04671587a975fb3c9fd2e0266336c96635e62cf6951  ma60/ma60d.f
ba0eb4e062db75e476f94bdc207d8f19e7bd6381f3a4107405e715b00a52ac05  mc71/mc71d.f
d073b0a641e48b4d7066c8a4bfff8bbbd065f41ac8b8d46691563cb1085a8ae7  fd15/fd15d.f
fe7806c3c9c7cebd32bfca09b4e8bdd86bd06da9045acf87a2017a184d15441c  hsl_mp48/hsl_mp48d.f90
2c063487f11a1122fb6048bf81bf2197b2f5e3599b5e88e870096c015ad60955  hsl_mc66/hsl_mc66d.f90
7a602665f29f81d5a0184e0261e67c25d926bcfba1f4329d1b4557c82cfbbdce  hsl_mc79/hsl_mc79i.f90
b4e4bf159ed82d1b70cb3ad8b36dfd2f31a6a80de046dcb278eb42e24371b4d4  hsl_mc79/C/hsl_mc79i_ciface.f90
bdb0e987bb2ce75c6ac4d8ab574098b766637a241a024eb5c0e2ff16ca1d7b53  hsl_mc79/C/hsl_mc79i.h
66667581aa88b99abc93d0f7151bc976faf7315b1d272ae97364d5bbae5c497f  hsl_mp01/hsl_mp01i.f90
7aa0efe88a8f382dbc079b97ea3e66cf6820bc459422f15e447a9a435f158f08  hsl_zd11/hsl_zd11d.f90
9531b8fb191f8436e7575baa0fdc60d48702e82d9fd4b874a342f7fbdeff64ce  hsl_fa14/hsl_fa14d.f90
0227347505491f61312c24129043a7fbc4a1bfc24550eb17b081b0b3ebb89bda  hsl_mc65/hsl_mc65d.f90
7cb52c0f07d08d8954257d3ce59d5a66f94976acffe42e109d062cd79af885a3  ma42/ma42d.f
a4318c5c5923b27984d847f1a7329518faee69a58d64d4641a6513c72c64701b  ma50/ma50d.f
b4240cf04cc85aa79563926d2c563fda2f7c346ebde8cdc821d114090b65fa01  ma52/ma52d.f
6f40e5872e9fcb62cba2fc49c246aeae13270dfe8159b9c78328f94b6ce7c151  mc13/mc13d.f
a46e2490f1760422579e0e4f74e90b1d058ce65346c3798eddd068f6b9b24224  mc21/mc21d.f
866badf468eabd0b8e57e746e1b56efd459fc73a72bd79310dd2cc0adfc9a931  mc29/mc29d.f
f40d712918a069435fd94a5c5df8ba579c0cfd8f9b91c821b6e26767d8e78aa6  mc46/mc46d.f
cc118acc23f6cebb89b1da4c18a741d4db97fd86ad43f91a0fde94ee7da67f07  mc59/mc59d.f"

if ! (cd "$ROOT" && echo "$MANIFEST" | sha256sum -c --quiet -); then
    if [ "${ALLOW_HSL_DRIFT:-0}" = "1" ]; then
        echo "stage_hsl.sh: WARNING: deck checksums differ from the verified" >&2
        echo "snapshot; proceeding because ALLOW_HSL_DRIFT=1." >&2
    else
        echo "stage_hsl.sh: libHSL $VERSION contains decks that differ from the" >&2
        echo "verified snapshot (2026.8.4). Diff the decks above against a" >&2
        echo "known-good build and re-run the parity kits before recording the" >&2
        echo "new checksums here. Set ALLOW_HSL_DRIFT=1 to build anyway." >&2
        exit 1
    fi
fi

stage() { mkdir -p "$(dirname "$2")" && cp "$ROOT/$1" "$2"; }

# ma48 / ma51: deck + dependency bundle -> shared libraries
stage ma48/ma48d.f "$STAGE/lib/ma48/src/ma48d.f"
cat "$ROOT"/mc59/mc59d.f "$ROOT"/ma50/ma50d.f "$ROOT"/mc13/mc13d.f \
    "$ROOT"/mc21/mc21d.f "$ROOT"/mc29/mc29d.f "$ROOT"/mc71/mc71d.f \
    > "$STAGE/lib/ma48/src/ddeps.f"
stage ma51/ma51d.f "$STAGE/lib/ma51/src/ma51d.f"
cat "$STAGE/lib/ma48/src/ddeps.f" "$ROOT"/ma48/ma48d.f \
    > "$STAGE/lib/ma51/src/ddeps.f"

# MP48: sources staged for the solver build (patched there); common90
# provides the hsl_mp01 module
stage hsl_mp48/hsl_mp48d.f90 "$STAGE/lib/mp48/src/hsl_mp48d.f90"
cat "$ROOT"/ma42/ma42d.f "$STAGE/lib/ma51/src/ddeps.f" \
    "$ROOT"/mc46/mc46d.f "$ROOT"/ma52/ma52d.f \
    > "$STAGE/lib/mp48/src/ddeps.f"
stage hsl_mp01/hsl_mp01i.f90 "$STAGE/lib/mp48/src/common90.f90"

# MC66 + its module dependencies
stage hsl_mc66/hsl_mc66d.f90 "$STAGE/lib/mc66/src/hsl_mc66d.f90"
cat "$ROOT"/hsl_zd11/hsl_zd11d.f90 "$ROOT"/hsl_fa14/hsl_fa14d.f90 \
    "$ROOT"/hsl_mc65/hsl_mc65d.f90 \
    > "$STAGE/lib/mc66/src/ddeps90.f90"

# MC79 (integer-only, dependency-free) + C interface
stage hsl_mc79/hsl_mc79i.f90          "$STAGE/lib/mc79/src/hsl_mc79i.f90"
stage hsl_mc79/C/hsl_mc79i_ciface.f90 "$STAGE/lib/mc79/src/hsl_mc79i_ciface.f90"
stage hsl_mc79/C/hsl_mc79i.h          "$STAGE/lib/mc79/include/hsl_mc79i.h"

# -condest diagnostics
stage ma60/ma60d.f "$STAGE/lib/ma60/src/ma60d.f"
stage mc71/mc71d.f "$STAGE/lib/mc71/src/mc71d.f"
stage fd15/fd15d.f "$STAGE/lib/fd15/fd15d.f"

# licence + provenance
mkdir -p "$STAGE/LICENSES/hsl"
cp "$ROOT/LICENCE" "$STAGE/LICENSES/hsl/LICENCE"
echo "libHSL $VERSION" > "$STAGE/LICENSES/hsl/VERSION"

echo "stage_hsl.sh: staged libHSL $VERSION into $STAGE"
