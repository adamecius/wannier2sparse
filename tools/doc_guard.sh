#!/usr/bin/env bash
# Pass iff the diff vs BASE changes only comments in C/C++ files (+ README/FINDINGS).
set -euo pipefail
BASE="${1:-origin/master}"
ALLOW_NONSRC='^(README\.md$|docs/FINDINGS\.md$|tools/doc_guard\.sh$)'
strip(){ gcc -fpreprocessed -dD -E -P -x c++ - 2>/dev/null; }  # remove comments, no expand
fail=0
while IFS= read -r f; do
  [ -z "$f" ] && continue
  case "$f" in
    *.hpp|*.cpp|*.h|*.cc|*.hh)
      if ! git cat-file -e "$BASE:$f" 2>/dev/null; then
        echo "GUARD: new source file '$f' not allowed"; fail=1; continue; fi
      if [ "$(git show "$BASE:$f" | strip)" != "$(git show "HEAD:$f" | strip)" ]; then
        echo "GUARD: non-comment change in $f"; fail=1; fi ;;
    *)
      echo "$f" | grep -qE "$ALLOW_NONSRC" || { echo "GUARD: disallowed file $f"; fail=1; } ;;
  esac
done < <(git diff --name-only "$BASE"...HEAD)
[ "$fail" -eq 0 ] && echo "GUARD: clean (comments + README only)." || { echo "GUARD: FAILED"; exit 1; }
