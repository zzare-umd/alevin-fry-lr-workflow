#!/usr/bin/env bash
# Builds the alevin-fry-lr binary from the long-read branch of alevin-fry into tool/bin/alevin-fry.
#
#   usage: tool/install.sh [install dir]      (default: the tool/ directory of this repository)
#
# Needs git, a Rust toolchain >= 1.85 (https://rustup.rs), cmake >= 3.5 (or CMAKE=/path/to/cmake) and
# a C compiler; cmake builds the compression library alevin-fry links. The build takes a few minutes.
# AF_REF pins the commit we tested; AF_REF=long-read builds the head of the branch.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${1:-$HERE}
AF_REF=${AF_REF:-d63db0c84066fa2a942fbaa603acaac9a29eb252}
for t in git cargo cc "${CMAKE:-cmake}"; do command -v "$t" > /dev/null || { echo "install.sh: '$t' not found in PATH (cmake can also be given as CMAKE=/path/to/cmake)" >&2; exit 1; }; done
if [ "$AF_REF" = long-read ]; then
  cargo install --locked --git https://github.com/COMBINE-lab/alevin-fry --branch long-read --root "$ROOT" alevin-fry
else
  cargo install --locked --git https://github.com/COMBINE-lab/alevin-fry --rev "$AF_REF" --root "$ROOT" alevin-fry
fi
echo "installed: $ROOT/bin/alevin-fry"
"$ROOT/bin/alevin-fry" --version
