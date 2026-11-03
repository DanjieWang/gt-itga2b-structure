#!/usr/bin/env Rscript
# =============================================================================
#  make_figure3.R
#  Structural modelling of two ITGA2B variants on the resting alphaIIb-beta3
#  ectodomain (PDB 3FCS) and rendering of Figure 3.
#
#  Inputs (expected in the working directory):
#    3FCS.pdb            - coordinates from RCSB PDB
#    panel_template.html - 3Dmol.js rendering template
#    3Dmol-min.js        - NOT included in this repository; download it first:
#                          curl -sL https://3Dmol.org/build/3Dmol-min.js -o 3Dmol-min.js
#
#  Outputs (in ./out):
#    A_Gln707_WT.png, B_Gln707_Ter.png, A_Gly267_WT.png, B_Gly267_skip.png
#    figure3_gln707_stop.png, figure3_exon7_skipping.png
#
#  Requirements: R with bio3d, jsonlite, magick; Microsoft Edge or Google Chrome
#                (used headlessly to screenshot the rendered panels).
#  Run:          Rscript make_figure3.R
# =============================================================================

suppressPackageStartupMessages({
  library(bio3d)
  library(jsonlite)
  library(magick)
})

BASE <- normalizePath(".", mustWork = TRUE)
setwd(BASE)
OUT <- file.path(BASE, "out")
dir.create(OUT, showWarnings = FALSE)
file.copy("3Dmol-min.js", file.path(OUT, "3Dmol-min.js"), overwrite = TRUE)

# ---- Locate a Chromium-based browser (Edge or Chrome) ------------------------
EDGE <- Sys.getenv("BROWSER_BIN", unset = NA)
if (is.na(EDGE)) {
  cand <- c("C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
            "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
            "C:/Program Files/Google/Chrome/Application/chrome.exe",
            "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/usr/bin/google-chrome", "/usr/bin/chromium")
  cand <- cand[file.exists(cand)]
  if (!length(cand)) {
    stop("No Chromium-based browser found. Set BROWSER_BIN to the full path of ",
         "msedge.exe or chrome.exe.")
  }
  EDGE <- cand[1]
}
cat("Browser:", EDGE, "\n")

# Viewer canvas size; the screenshot is taken at exactly this size.
VW <- 1400
VH <- 1180

# ---------------------------------------------------------------- 1. Structure
pdb <- read.pdb("3FCS.pdb", verbose = FALSE)
at  <- pdb$atom
pdb_text <- paste(readLines("3FCS.pdb", warn = FALSE), collapse = "\n")

AA1 <- c(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H",
         ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W",
         TYR="Y", VAL="V")
AA3 <- c(ALA="Ala", ARG="Arg", ASN="Asn", ASP="Asp", CYS="Cys", GLN="Gln", GLU="Glu",
         GLY="Gly", HIS="His", ILE="Ile", LEU="Leu", LYS="Lys", MET="Met", PHE="Phe",
         PRO="Pro", SER="Ser", THR="Thr", TRP="Trp", TYR="Tyr", VAL="Val")

# Numbering check: 3FCS chain A residue numbers are 31 lower than NP_000410.
seqres <- grep("^SEQRES", readLines("3FCS.pdb", warn = FALSE), value = TRUE)
chain_a <- unlist(lapply(seqres[substr(seqres, 12, 12) == "A"], function(l) {
  r <- trimws(substring(l, seq(20, 68, 4), seq(22, 70, 4)))
  r[r != ""]
}))
chain_a <- paste(unname(AA1[chain_a]), collapse = "")
stopifnot(substr(chain_a, 676, 676) == "Q", substr(chain_a, 236, 236) == "G")
cat("Numbering check passed: 676 = Q, 236 = G\n")

# ------------------------------------------------------------ 2. Local context
# Residues whose CA lies within 13 A of the site: drawn as a cartoon ribbon
# coloured by secondary structure, to give the panel some structural context.
cartoon_res <- function(at, chain, resno, radius = 13) {
  ca <- at[at$chain == chain & at$elety == "CA", ]
  c0 <- ca[ca$resno == resno, ]
  d  <- sqrt((ca$x - c0$x)^2 + (ca$y - c0$y)^2 + (ca$z - c0$z)^2)
  sort(unique(ca$resno[d <= radius & ca$resno != resno]))
}

# Neighbouring residues (any atom within `cutoff`) and hydrogen bonds.
# `deleted`: residues removed first, to mimic a truncation or an exon skip.
env_panel <- function(at, chain, resno, cutoff = 4.5, hb_max = 3.5, deleted = integer(0)) {
  if (length(deleted)) at <- at[!(at$chain == "A" & at$resno %in% deleted), ]
  center <- at[at$chain == chain & at$resno == resno, ]
  other  <- at[!(at$chain == chain & at$resno == resno), ]

  # Distance from each atom of the residue of interest to every other atom.
  # Do not use dist() on 17,000 atoms: that matrix needs ~2 GB.
  sub <- sqrt(outer(center$x, other$x, "-")^2 +
              outer(center$y, other$y, "-")^2 +
              outer(center$z, other$z, "-")^2)

  hit  <- which(sub <= cutoff, arr.ind = TRUE)
  near <- unique(other[hit[, "col"], c("chain", "resno", "resid")])
  near <- near[order(near$chain, near$resno), ]

  # Hydrogen bond: both atoms are N or O, 2.0 < d <= 3.5 A, and the pair is not
  # a peptide bond (same chain, adjacent residues).
  cno <- which(substr(center$elety, 1, 1) %in% c("N", "O"))
  ono <- which(substr(other$elety,  1, 1) %in% c("N", "O"))
  hb <- list()
  for (i in cno) for (j in ono) {
    d <- sub[i, j]
    if (d <= 2.0 || d > hb_max) next
    if (center$chain[i] == other$chain[j] && abs(center$resno[i] - other$resno[j]) == 1) next
    hb[[length(hb) + 1]] <- list(
      x1 = round(center$x[i], 3), y1 = round(center$y[i], 3), z1 = round(center$z[i], 3),
      x2 = round(other$x[j], 3),  y2 = round(other$y[j], 3),  z2 = round(other$z[j], 3),
      d  = round(d, 2))
  }
  list(near = near, hb = hb)
}

# Residue data frame -> [[chain, resno, resname], ...] for JSON.
to_list <- function(df) lapply(seq_len(nrow(df)), function(i)
  list(df$chain[i], as.integer(df$resno[i]), df$resid[i]))

# Nearest residue to a coordinate, used to label both ends of a hydrogen bond.
nearest_res <- function(at, x, y, z) {
  d <- (at$x - x)^2 + (at$y - y)^2 + (at$z - z)^2
  at[which.min(d), c("chain", "resno", "resid")]
}

# 2.1 Wild-type Gln707 (3FCS residue 676)
wt707 <- env_panel(at, "A", 676)

# 2.2 p.Gln707*: keep residues 1-706 only
trunc_del <- 676:max(at$resno[at$chain == "A"])
trunc     <- env_panel(at, "A", 675, deleted = trunc_del)

# 2.3 Wild-type Gly267 (3FCS residue 236, last codon of exon 7)
wt267 <- env_panel(at, "A", 236)

# 2.4 Exon 7 skipping: delete 3FCS 193-235 (= NP_000410 224-266, 43 residues)
skip_del <- 193:235
skip     <- env_panel(at, "A", 236, deleted = skip_del)

# 2.5 New backbone junction: the two CA atoms flanking the deletion
sel_skip <- at[!(at$chain == "A" & at$resno %in% skip_del), ]
ca1 <- sel_skip[sel_skip$chain == "A" & sel_skip$resno == 192 & sel_skip$elety == "CA", ]
ca2 <- sel_skip[sel_skip$chain == "A" & sel_skip$resno == 236 & sel_skip$elety == "CA", ]
junction <- list(x1 = round(ca1$x[1], 3), y1 = round(ca1$y[1], 3), z1 = round(ca1$z[1], 3),
                 x2 = round(ca2$x[1], 3), y2 = round(ca2$y[1], 3), z2 = round(ca2$z[1], 3))

# 2.6 Residues to label: the partners at the far end of each hydrogen bond
hblab <- function(cfg, center) {
  labs <- lapply(cfg$hb, function(h) {
    rbind(nearest_res(at, h$x1, h$y1, h$z1), nearest_res(at, h$x2, h$y2, h$z2))
  })
  labs <- do.call(rbind, labs)
  labs <- unique(labs[labs$resno != center, c("chain", "resno")])
  lapply(seq_len(nrow(labs)), function(i) list(labs$chain[i], as.integer(labs$resno[i])))
}

# 2.7 Panel configurations
mk_cfg <- function(center, near, hb, labels, mainlabel, deleted = integer(0),
                   ghost = NULL, junction = NULL, dx = 2, dy_label = 1.5) {
  cfg <- list(center = center, deleted = as.integer(deleted), near = to_list(near),
              hb = hb, labels = labels, mainlabel = mainlabel, dx = dx, dy_label = dy_label)
  if (!is.null(ghost))    cfg$ghost <- ghost
  if (!is.null(junction)) cfg$junction <- junction
  cfg
}

lab707 <- hblab(wt707, 676)
lab267 <- hblab(wt267, 236)
# Tyr265 is itself encoded by exon 7, so it cannot be labelled in the skip model
lab267b <- Filter(function(r) !(r[[2]] %in% skip_del), lab267)

# Panel a uses the wild-type neighbourhood; panel b reuses the same `near` list
# so that both panels share one camera and the two structures stay aligned.
cfgA707 <- mk_cfg(676, wt707$near, wt707$hb, lab707, "Gln707")

# Truncation panel: centred on the new C-terminus Asn706 (3FCS 675); deleted
# neighbours are drawn as pale grey ghosts.
cfgB707 <- mk_cfg(675, wt707$near, list(), list(), "Asn706",
                  trunc_del,
                  ghost = to_list(wt707$near[wt707$near$resno %in% trunc_del, ]),
                  dx = 4.5, dy_label = -2)

cfgA267 <- mk_cfg(236, wt267$near, wt267$hb, lab267, "Gly267")
cfgB267 <- mk_cfg(236, wt267$near, skip$hb, lab267b, "Gly267", skip_del,
                  ghost = to_list(wt267$near[wt267$near$resno %in% skip_del, ]),
                  junction = junction)

# The 2.87 A label overlaps the 3.01 A label on screen; nudge only the label
shift_287 <- function(cfg) {
  cfg$hb <- lapply(cfg$hb, function(h) { if (h$d == 2.87) h$dy <- 1.4; h })
  cfg
}

# Both panels are framed on the cartoon ribbon plus the residue of interest,
# so the field of view is identical in a and b.
add_cartoon <- function(cfg, resno) {
  rs <- cartoon_res(at, "A", resno)
  cfg$cartoon <- lapply(rs, function(r) list("A", as.integer(r)))
  cfg$fit <- c(cfg$cartoon, list(list("A", as.integer(resno))))
  cfg
}
cfgA707 <- add_cartoon(cfgA707, 676)
cfgB707$cartoon <- cfgA707$cartoon
cfgB707$fit     <- cfgA707$fit
cfgA267 <- add_cartoon(cfgA267, 236)
cfgB267$cartoon <- cfgA267$cartoon
cfgB267$fit     <- cfgA267$fit

cfgA267 <- shift_287(cfgA267)
cfgB267$fit <- cfgA267$fit
cfgB267 <- shift_287(cfgB267)

# ---------------------------------------------------------------- 3. Write HTML
tpl <- paste(readLines("panel_template.html", warn = FALSE), collapse = "\n")
aa3_js <- toJSON(as.list(AA3), auto_unbox = TRUE)

build_html <- function(cfg, name) {
  js <- toJSON(cfg, auto_unbox = TRUE, digits = 3, null = "null")
  h  <- sub("{{PDB}}", pdb_text, tpl, fixed = TRUE)
  h  <- sub("{{AA3}}", aa3_js,  h,   fixed = TRUE)
  h  <- sub("{{CFG}}", js,      h,   fixed = TRUE)
  h  <- sub("{{W}}",   VW,      h,   fixed = TRUE)
  h  <- sub("{{H}}",   VH,      h,   fixed = TRUE)
  f  <- file.path(OUT, paste0(name, ".html"))
  writeLines(h, f, useBytes = TRUE)
  f
}

# ---------------------------------------------------------------- 4. Screenshot
# A fresh --user-data-dir per shot; otherwise Edge serves the cached page and
# the screenshot does not reflect the current configuration.
shot <- function(html, png) {
  png  <- normalizePath(png, mustWork = FALSE)
  prof <- file.path(tempdir(), paste0("edge-", as.integer(runif(1, 1, 1e9))))
  url  <- paste0("file:///", gsub("\\\\", "/", normalizePath(html)))
  args <- c("--headless=new", "--disable-gpu", "--hide-scrollbars",
            "--force-device-scale-factor=1",
            sprintf('--user-data-dir="%s"', prof),
            sprintf('--screenshot="%s"', png),
            sprintf("--window-size=%d,%d", VW, VH),
            "--virtual-time-budget=20000",
            sprintf('"%s"', url))
  system(paste(shQuote(EDGE), paste(args, collapse = " ")), intern = FALSE, wait = TRUE)
  if (!file.exists(png) || file.size(png) < 5000) {
    stop("Screenshot failed: ", png,
         "\nCheck the browser path and that 3Dmol-min.js is in out/.")
  }
  png
}

panels <- list(A_Gln707_WT = cfgA707, B_Gln707_Ter = cfgB707,
               A_Gly267_WT = cfgA267, B_Gly267_skip = cfgB267)
for (nm in names(panels)) {
  shot(build_html(panels[[nm]], nm), file.path(OUT, paste0(nm, ".png")))
  cat("rendered:", nm, "\n")
}

# ------------------------------------------------------ 5. Composite the figure
PW <- 820; PH <- 680          # panel size inside the composite
CW <- 1700; CH <- 960         # composite canvas

# The two panels are not cropped individually: they share one camera, so scaling
# the whole canvas keeps the molecules in register.
prep_pair <- function(fa, fb) {
  fit <- function(f) {
    im <- image_read(file.path(OUT, paste0(f, ".png")))
    im <- image_resize(im, sprintf("%dx%d", PW, PH))
    image_extent(im, sprintf("%dx%d", PW, PH), gravity = "center", color = "white")
  }
  list(fit(fa), fit(fb))
}

txt <- function(canvas, s, size, bold = FALSE, x, y, color = "#3C3C3C") {
  image_annotate(canvas, s, font = "Arial", size = size,
                 weight = if (bold) 700 else 400, color = color,
                 location = sprintf("+%d+%d", x, y), gravity = "northwest")
}

compose <- function(out, title, subtitle, cap, pa, pb) {
  canvas <- image_blank(CW, CH, "white")
  canvas <- image_composite(canvas, pa, offset = "+20+120")
  canvas <- image_composite(canvas, pb, offset = "+860+120")
  canvas <- txt(canvas, title,    30, TRUE,  20, 14, "#141414")
  canvas <- txt(canvas, subtitle, 22, FALSE, 20, 56, "#3C3C3C")
  canvas <- txt(canvas, "a", 34, TRUE, 26, 86, "#141414")
  canvas <- txt(canvas, "b", 34, TRUE, 866, 86, "#141414")
  y <- 812
  for (ln in strwrap(cap, width = 138)) {
    canvas <- txt(canvas, ln, 19, FALSE, 20, y, "#464646")
    y <- y + 24
  }
  junk <- image_write(canvas, file.path(OUT, out))
  cat("saved:", out, "\n")
  canvas
}

p707 <- prep_pair("A_Gln707_WT", "B_Gln707_Ter")
p267 <- prep_pair("A_Gly267_WT",  "B_Gly267_skip")

# Also keep the aligned single panels, which are handy for slide decks.
junk <- image_write(p707[[1]], file.path(OUT, "panel_A_Gln707_WT.png"))
junk <- image_write(p707[[2]], file.path(OUT, "panel_B_Gln707_Ter.png"))
junk <- image_write(p267[[1]], file.path(OUT, "panel_A_Gly267_WT.png"))
junk <- image_write(p267[[2]], file.path(OUT, "panel_B_Gly267_skip.png"))

cap707 <- paste("a, Wild-type alphaIIb around Gln707 (green). A single predicted hydrogen bond of 3.24 A links its side-chain oxygen to Arg714.",
                "b, c.2119C>T (p.Gln707*) introduces a premature stop codon, so the model ends at Asn706 (green, b); all residues from 707 onward,",
                "including the remainder of calf-1, calf-2 and the transmembrane/cytoplasmic domains, are lost. Residues are numbered as in NP_000410;",
                "the corresponding PDB 3FCS numbers are 31 lower (Gln707 = residue 676).")
cap267 <- paste("a, Wild-type alphaIIb around Gly267 (green). Its exposed backbone N and O atoms form four predicted hydrogen bonds (2.87-3.15 A)",
                "with Ser203, Tyr265, Ser269 and Gly287.",
                "b, c.799+1G>A abolishes the exon 7 splice donor and is predicted to cause skipping of exon 7 (43 residues, 224-266); the residues encoded by",
                "exon 7 are drawn in pale grey and the new backbone junction is marked by the purple dashed line.")

compose("figure3_gln707_stop.png",
        "ITGA2B c.2119C>T (p.Gln707*), homozygous",
        "Wild type (a) vs. predicted truncation (b)", cap707, p707[[1]], p707[[2]])

compose("figure3_exon7_skipping.png",
        "ITGA2B c.799+1G>A, homozygous",
        "Wild type (a) vs. predicted exon 7 skipping (b)", cap267, p267[[1]], p267[[2]])

cat("\nDone. Four single panels and two composites are in:", OUT, "\n")
