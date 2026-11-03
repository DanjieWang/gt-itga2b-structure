# Structural modelling of two ITGA2B variants (p.Gln707\* and c.799+1G>A)

Code that reproduces Figure 3 of the study *"Two novel homozygous ITGA2B
variants associated with type I Glanzmann thrombasthenia"* (manuscript under
review).

Two ITGA2B variants are modelled on the crystal structure of the resting
αIIbβ3 ectodomain (PDB [3FCS](https://www.rcsb.org/structure/3FCS)):

- **c.2119C>T (p.Gln707\*)** — premature termination; residues 707–1039 are
  removed from the model and the new C-terminus (Asn706) is highlighted.
- **c.799+1G>A** — predicted in-frame skipping of exon 7; residues 224–266
  (NP_000410 numbering) are removed and the new backbone junction is drawn as
  a purple dashed line.

For each site the local environment is extracted (any atom within 4.5 Å) and
hydrogen bonds are assigned by a geometric criterion (N/O···N/O,
2.0–3.5 Å, excluding peptide-bonded neighbours). Panels are rendered with
[3Dmol.js](https://3dmol.org), screenshotted headlessly and composited with
ImageMagick.

## Contents

| File | Description |
|---|---|
| `make_figure3.R` | The whole analysis and figure pipeline |
| `panel_template.html` | Rendering template (3Dmol.js scene definition) |
| `3FCS.pdb` | Input coordinates, downloaded from RCSB PDB |
| `LICENSE` | MIT licence |

## Requirements

- R ≥ 4.2 with **bio3d**, **jsonlite**, **magick**
  `install.packages(c("bio3d", "jsonlite", "magick"), type = "binary")`
- Microsoft Edge or Google Chrome, used in headless mode for screenshots
  (auto-detected; override with the environment variable `BROWSER_BIN`)
- **3Dmol.js** is *not* included here. Fetch it next to the script:
  `curl -sL https://3Dmol.org/build/3Dmol-min.js -o 3Dmol-min.js`

Keep the checkout in a path without non-ASCII characters: on a Chinese-locale
Windows install R cannot translate such paths for the browser and ImageMagick
calls.

## Usage

```bash
Rscript make_figure3.R
```

Output goes to `out/`: four single panels (`A_Gln707_WT.png`,
`B_Gln707_Ter.png`, `A_Gly267_WT.png`, `B_Gly267_skip.png`), the aligned
versions of the same panels (`panel_*.png`) and the two composites
(`figure3_gln707_stop.png`, `figure3_exon7_skipping.png`). Rendering takes
about one minute.

## Methodological notes

1. **Residue numbering** — the text follows NP_000410 (ITGA2B); 3FCS chain A
   numbering is 31 lower (Gln707 → residue 676, Gly267 → residue 236). The
   offset is verified against the SEQRES records at start-up, and the script
   stops if the check fails.
2. **Hydrogen bonds** — assigned from heavy-atom geometry only (N/O pairs,
   2.0–3.5 Å, peptide-bonded neighbours excluded). Donor–acceptor angles and
   proton positions are not considered, so these are predicted, not
   experimentally determined, bonds.
3. **Mutant models** — the truncation and the exon 7 skip are schematic
   predictions obtained by deleting the affected residues from the wild-type
   structure; they are not experimental structures. The deleted residues are
   drawn in pale grey at their wild-type positions for orientation.
4. **Reproducibility** — panels a and b share one camera and one colour scale,
   so the two are directly comparable.

## Licence

Copyright (c) 2026 Danjie Wang. Released under the [MIT Licence](LICENSE).
