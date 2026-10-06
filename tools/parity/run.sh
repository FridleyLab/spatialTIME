#!/usr/bin/env bash
#
# Regenerate the cross-branch parity fixture.
#
#   bash tools/parity/run.sh
#
# Captures spatialTIME 1.4.0 (the feature-alex tip) and the current working tree in
# two separate R processes, diffs them against the contract in compare.R, and writes
# tests/testthat/fixtures/parity-v1.4.0.rds.
#
# Exits non-zero if any measured result contradicts its declared verdict, so this is
# a check and not only a generator. Run it after any change to the statistics, and
# read the verdict table rather than only the exit code -- a column moving from
# MUST_MATCH to MUST_DIFFER is a decision, not a detail.
#
# Nothing here is part of the built package: `^tools$` is in .Rbuildignore, and the
# fixture is too, so `R CMD check` on the tarball never sees either and
# test-parity-v1.4.0.R skips. From a git checkout `devtools::test()` runs it.

set -euo pipefail

# The 1.4.0 side. A commit, not a branch name, so this keeps working if feature-alex
# moves or is deleted.
OLD_REF="${OLD_REF:-20960eb}"

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

WT="$(mktemp -d /tmp/st-alex.XXXXXX)"
OUT="$(mktemp -d /tmp/st-parity.XXXXXX)"
cleanup() {
  git worktree remove --force "$WT" 2>/dev/null || true
  rm -rf "$OUT"
}
trap cleanup EXIT

echo "==> worktree for $OLD_REF at $WT"
git worktree add --detach "$WT" "$OLD_REF" >/dev/null

# --pkg="$REPO" is deliberately the LIVE working tree, not a copy. For a pure-R
# package pkgload::load_all() writes nothing into the package directory, and a copy
# would only add a way for the two to drift. Note that git plumbing cannot snapshot
# this side anyway: `git stash create` has no --include-untracked, and the
# disk-backed store files were untracked when this harness was written.
echo "==> capturing old (1.4.0)"
Rscript tools/parity/capture.R --pkg="$WT"   --side=old --out="$OUT/old.rds"
echo "==> capturing new (working tree)"
Rscript tools/parity/capture.R --pkg="$REPO" --side=new --out="$OUT/new.rds"

echo "==> comparing"
Rscript tools/parity/compare.R "$OUT/old.rds" "$OUT/new.rds" \
  --fixture=tests/testthat/fixtures/parity-v1.4.0.rds

echo
echo "==> done. Review the verdict table above, then:"
echo "    Rscript -e 'devtools::test(filter = \"parity\")'"
