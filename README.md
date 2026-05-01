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

```bash
python3 freebayes_somatic_subset.py --input-dir <raw-vcf-dir> -o docs/default \
    --min-tumor-vaf 0.05 --max-normal-vaf 0.02 --qual 20  # default cutoffs
```

After regenerating, refresh `docs/default/manifest.json` so the explorer
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
