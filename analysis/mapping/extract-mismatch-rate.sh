#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob  # ensures empty glob expands to nothing, not literal *.cram.stats

OUT="mismatches.tsv"
echo -e "sample_id\tmismatches" > "$OUT"

for dir in */ ; do
  [[ -d "$dir" ]] || continue
  sample_id="${dir%/}"

  # Pick up samtools stats files (assuming only one per sample)
  stats_files=("$dir"/*.cram.stats)

  if [[ ${#stats_files[@]} -gt 0 ]]; then
    stats_file="${stats_files[0]}"
    mismatches="$(awk -F $'\t' '$1=="SN" && $2=="error rate:" {print $3}' "$stats_file" | head -n1)"
    [[ -n "$mismatches" ]] || mismatches="NA"
  else
    mismatches="NA"
  fi

  printf "%s\t%s\n" "$sample_id" "$mismatches" >> "$OUT"
done

echo "✅ Wrote $(($(wc -l < "$OUT") - 1)) sample entries to $OUT"
