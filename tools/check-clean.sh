#!/bin/sh

set -eu

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
ROOT_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"

cd "$ROOT_DIR"

./cleanup

R CMD build . >/tmp/simplicialdepth_build.log 2>&1
TARBALL=$(awk '/^\* building/{print $NF}' /tmp/simplicialdepth_build.log | tail -n1 | tr -d '‘’')

if [ -z "$TARBALL" ] || [ ! -f "$TARBALL" ]; then
  echo "Failed to locate built tarball from R CMD build output." 1>&2
  exit 1
fi

# Default mode is --as-cran, but callers can pass custom check args.
if [ "$#" -eq 0 ]; then
  CHECK_ARGS="--as-cran"
else
  CHECK_ARGS="$*"
fi

# shellcheck disable=SC2086
R CMD check $CHECK_ARGS "$TARBALL"
