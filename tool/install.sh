#!/usr/bin/env bash
# Builds the alevin-fry-lr binary: the long-read branch of alevin-fry together with the libradicl
# it depends on (a path dependency, so both are cloned side by side and libradicl gets a small patch).
#
#   usage: tool/install.sh [build dir]      (default: the tool/ directory of this repository)
#   result: <build dir>/alevin-fry/target/release/alevin-fry
#
# Needs: git, a Rust toolchain >= 1.85 (https://rustup.rs), cmake >= 3.5 and a C compiler
# (cmake builds the zlib-ng that alevin-fry links). Takes a few minutes; the build uses fat LTO.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PREFIX=${1:-$HERE}
AF_REF=${AF_REF:-a9c9b4265342d1259add1567f83d0eeda254e3a5}   # alevin-fry, branch long-read
LR_REF=${LR_REF:-6cd08c35085be81783042694dab4bce33f828def}   # libradicl, branch develop (v0.13.0)
mkdir -p "$PREFIX"; cd "$PREFIX"
for t in git cargo cc "${CMAKE:-cmake}"; do command -v "$t" > /dev/null || { echo "install.sh: '$t' not found in PATH (cmake can also be given as CMAKE=/path/to/cmake)" >&2; exit 1; }; done
if [ ! -d libradicl ]; then git clone --quiet https://github.com/COMBINE-lab/libradicl.git; fi
( cd libradicl && git fetch --quiet origin && git checkout --quiet "$LR_REF" \
  && if git apply --check "$HERE/libradicl-long-read.patch" 2> /dev/null; then git apply "$HERE/libradicl-long-read.patch"; echo "libradicl: patch applied";
     elif grep -q skip_read_name libradicl/src/record.rs; then echo "libradicl: patch already applied"; else echo "libradicl: the patch does not apply to $LR_REF" >&2; exit 1; fi )
if [ ! -d alevin-fry ]; then git clone --quiet --branch long-read https://github.com/COMBINE-lab/alevin-fry.git; fi
( cd alevin-fry && git fetch --quiet origin long-read && git checkout --quiet "$AF_REF" )
( cd alevin-fry && cargo build --release )
echo "built: $PREFIX/alevin-fry/target/release/alevin-fry"
"$PREFIX/alevin-fry/target/release/alevin-fry" --version
