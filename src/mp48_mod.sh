#!/bin/bash
# Apply TEEMS modifications to the HSL MP48 sources staged in this
# directory (hsl_mp48d.f90, ddeps.f):
#   - 64-bit duplicate-detection work arrays (MYIW1/MYLMX) and CGLOB
#     allocation fix in hsl_mp48d.f90
#   - MA48/MA50/MC13/MC21/MC29/MC59/MC71 entry points renamed to Z*
#     so the statically patched copies cannot clash with the symbols
#     in the dynamically linked libma48/libma51
# The patches were generated from HSL MP48 2.1.1; `patch` fails loudly
# if the upstream sources ever drift.
set -e
cd "$(dirname "$0")"
patch --forward hsl_mp48d.f90 patches/hsl_mp48d.patch
patch --forward ddeps.f patches/ddeps.patch
