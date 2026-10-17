#!/bin/bash
# Track 3 specification v2 (English translation of the original prompts of the first two tracks, with the per-branch guidance removed): one codex exec per PR (gpt-6-astra / medium, may read GitHub online), 3 in parallel; skip PRs that already have a spec.md
B="$(cd "$(dirname "$0")" && pwd)"
one() { d="$1"; echo "start $d $(date +%T)"
  codex exec -m gpt-6-astra -c model_reasoning_effort=medium -s workspace-write -c sandbox_workspace_write.network_access=true --skip-git-repo-check -C "$d" "$(cat "$B/prompt.txt")" < /dev/null > "$d/codex.log" 2>&1
  echo "end $d rc=$? $([ -s "$d/spec.md" ] && echo ok || echo MISSING) $(date +%T)"; }
for d in "$B"/*/; do [ -s "$d/spec.md" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 5; done; one "$d" & done; wait
echo ALL_DONE
