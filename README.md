# ucn1

FreeBayes tumor-vs-normal somatic-variant exploration. Companion VCFs are
filtered into a high-confidence somatic-only subset by `freebayes_somatic_subset.py`,
and the static `docs/freebayes_explorer.html` browser is the visual front-end:
a single-page tool that loads multiple subset VCFs, supports per-site filters,
BED-region restriction, pairwise comparisons (Venn + counts), IGV linking, and
TSV export.

The site is hosted on GitHub Pages from the `docs/` directory.

## Live site

After GitHub Pages is enabled, the explorer is at:

> **https://walaj.github.io/ucn1/**

When the page loads it auto-fetches every VCF listed in `docs/default/manifest.json`
so the viewer comes up populated. Override the source directory by appending
`?dir=some/other/folder` to the URL.

## Repository layout

```
ucn1/
├── README.md
├── .gitignore
├── freebayes_somatic_subset.py        # Python script that produces the subset VCFs
└── docs/                              # GitHub Pages source
    ├── index.html                     # redirect → freebayes_explorer.html
    ├── freebayes_explorer.html        # the explorer (single-file app)
    └── default/                       # subset VCFs served by Pages
        ├── manifest.json              # list of VCFs the explorer auto-loads
        ├── *.somatic_subset.vcf
        └── somatic_subset_summary.tsv
```

## Running locally

Open `docs/freebayes_explorer.html` directly in a browser (drag VCFs onto the
upload box), **or** serve the `docs/` folder over HTTP so auto-load works:

```bash
python3 -m http.server 8000 --directory docs
# then visit http://localhost:8000/
```

Either path supports loading additional VCFs by drag-and-drop or via the
folder picker.

## Regenerating the subsets

The Python subsetter has three modes via `--preset`:

| Preset       | Intent                                  | Key cutoffs                                                                        | Output suffix              |
|--------------|-----------------------------------------|------------------------------------------------------------------------------------|----------------------------|
| `strict`     | High-specificity nomination set         | qual≥200, mqm≥60, T-VAF≥0.15, T-AD≥10, N-VAF≤0.005, N-AD=0, SNVs only              | `.strict_subset.vcf`       |
| `sensitive`  | Broad evidence net for rescue checks    | qual≥1, mqm≥30, T-VAF≥0.02, T-AD≥2, normal unrestricted, all variant types         | `.sensitive_subset.vcf`    |
| `custom` (default) | Use individual `--min-*/--max-*` flags as set on the CLI | (whatever you pass)                                                          | `.somatic_subset.vcf`      |

`strict ⊂ sensitive` by construction (looser cutoffs admit a superset). Each
output VCF gets a `##sarek_subset_preset=…` header line so the explorer can
auto-pair strict/sensitive partners.

### Single-preset run (default cutoffs)

```bash
python3 freebayes_somatic_subset.py --input-dir <raw-vcf-dir> -o docs/default \
    --min-tumor-vaf 0.05 --max-normal-vaf 0.02 --qual 20
```

### Strict + sensitive pair (recommended for rescue mode)

For each raw FreeBayes T/N VCF, produce both presets so the explorer can do
asymmetric comparisons (nominate from A.strict, rescue against B.sensitive
and vice-versa):

```bash
for vcf in *.freebayes.vcf; do
  python3 freebayes_somatic_subset.py "$vcf" -o docs/default --preset strict
  python3 freebayes_somatic_subset.py "$vcf" -o docs/default --preset sensitive
done
```

Each input produces two outputs in `docs/default/`:
`<stem>.strict_subset.vcf` and `<stem>.sensitive_subset.vcf`.

### Refresh the manifest

After regenerating, rebuild `docs/default/manifest.json` so the explorer
auto-loads the new files:

```bash
python3 - <<'PY'
import json, os
d = "docs/default"
files = sorted(f for f in os.listdir(d) if f.endswith(".vcf") or f.endswith(".vcf.gz"))
json.dump({"name":"ucn1 default subsets",
          "description":"FreeBayes tumor-vs-normal somatic subsets.",
          "files": files},
          open(os.path.join(d, "manifest.json"), "w"), indent=2)
PY
```

## Rescue mode (asymmetric comparison)

When comparing two samples A and B, the obvious-but-wrong question is "which
strict calls in A also pass strict in B?" That penalises real-but-marginal
variants in B that are buried under the strict cutoffs. The right question is:

  *"Which calls strictly nominated in A have **at least some evidence** in B?"*

This is the rescue: keep the high-specificity set as the nomination list, but
check membership against a high-sensitivity set on the partner side. Set
operations become asymmetric:

  - `A only` = `A.strict − B.sensitive` (in A confidently, **no evidence at all** in B)
  - `B only` = `B.strict − A.sensitive`
  - `A ∩ B` = `(A.strict ∩ B.sensitive) ∪ (B.strict ∩ A.sensitive)` (confident in one side, evidence in the other)
  - `A ∪ B` = `A.strict ∪ B.strict`

In the explorer, click the "rescue mode" chip in the comparison panel. Two
extra dropdowns appear for the sensitive partners; the **auto-pair** button
matches each side with the same-tumor-sample file of the opposite preset
(detected from the `##sarek_subset_preset=…` header), so if you've loaded a
strict + sensitive pair per sample the partners populate themselves.

## Enabling GitHub Pages (one-time)

1. Push this repo to GitHub: `git push -u origin main`.
2. On GitHub, go to **Settings → Pages**.
3. Under **Build and deployment** set **Source** to **Deploy from a branch**,
   **Branch** to `main`, **Folder** to `/docs`.
4. Wait ~30 seconds for the first build, then load
   `https://walaj.github.io/ucn1/`.

GitHub Pages serves `/docs/index.html` at the root URL; that page redirects
to `freebayes_explorer.html`, which in turn fetches `default/manifest.json`
and pulls each listed VCF.
