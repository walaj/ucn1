# 1. Generate strict + sensitive for each raw VCF (handle the 12 GB bladder_n1 last;
#    each preset run on it is the slow one)
for vcf in ~/Dropbox/AndyS/sarek/*.vcf; do
  python3 freebayes_somatic_subset.py "$vcf" -o docs/default --preset strict
  python3 freebayes_somatic_subset.py "$vcf" -o docs/default --preset sensitive
done

# 2. Refresh manifest.json
python3 - <<'PY'
import json, os
d = "docs/default"
files = sorted(f for f in os.listdir(d) if f.endswith(".vcf") or f.endswith(".vcf.gz"))
json.dump({"name":"ucn1 default subsets",
          "description":"FreeBayes tumor-vs-normal somatic subsets.",
          "files": files},
          open(os.path.join(d, "manifest.json"), "w"), indent=2)
PY

