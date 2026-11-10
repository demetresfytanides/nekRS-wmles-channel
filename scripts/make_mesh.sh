#!/bin/sh
# Build channel.re2 (and channel.ma2) from case/channel.box with Nek5000's
# genbox/genmap. Run from the directory holding the case files.
#   NEK_TOOLS=/path/to/Nek5000/bin  ../scripts/make_mesh.sh
set -e
CASE=channel
BIN="${NEK_TOOLS:+$NEK_TOOLS/}"
printf '%s.box\n' "$CASE" | "${BIN}genbox"
mv -f box.re2 "$CASE.re2"
printf '%s\n0.2\n' "$CASE" | "${BIN}genmap"
