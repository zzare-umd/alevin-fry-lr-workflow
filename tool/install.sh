#!/usr/bin/env bash
# Installs the long-read alevin-fry into tool/bin/alevin-fry (the path the workflow config uses).
# Needs a Rust toolchain >= 1.85 (https://rustup.rs), cmake and a C compiler; takes a few minutes.
set -euo pipefail
cd "$(dirname "$0")"
cargo install --locked --git https://github.com/COMBINE-lab/alevin-fry --branch long-read --root . alevin-fry
./bin/alevin-fry --version
