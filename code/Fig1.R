# ---------------------------------------------------------------------------
# 20260604_JS_AmylEx_Fig1.R — Amyloid Explorer database composition figure
#   (rebuild of the original Fig1 against the redesigned pipeline)
#
# Reads paper_figures/data/amyloid_explorer_db.tsv (the public 19-col view
# of the database) and produces a 7-panel figure matching the original
# Fig1 layout, updated for the redesigned pipeline:
#
#   A — 4 info cards (Structures, Proteins, Species, Techniques)
#   B — donut chart of top proteins by structure count
#   C — cryo-EM resolution histogram + median
#   D — Role pie (Pathological / Functional / Synthetic — the new S
#       category gets its own slice; old version was just P/F)
#   E — species bars (Human / Mouse / Other; normalised from the
#       UniProt-style long organism names)
#   F — sample-origin bars (Patient / Animal / Seeded / Recombinant —
#       all four sample_origin categories; old version dropped
#       Recombinant entirely)
#   G — technique bars (cryoEM / solidNMR / solutNMR; normalised from
#       "ELECTRON MICROSCOPY" etc.)
#
# Output: paper_figures/output/Fig1_v2.pdf
#
# IMPORTANT — cohort scoping:
#   This figure describes the SUBSET of the Amyloid Explorer archive
#   that went through the energetics analysis (Figs 2 onward) — i.e.,
#   the 775 PDBs that survived Nikos's RepairPDB + chain-fit QC.  Every
#   panel here is filtered to those PDBs.  The full public archive on
#   amyloidexplorer.org is larger (829 at the time of this rebuild);
#   the 54 archive-only structures are documented in the supplementary
#   exclusion table.  Reasoning: every figure in the paper should
#   reflect the same cohort so reviewers don't have to mentally bridge
#   numbers across figures.
#
# Differences from the old version (paper_figures/old/Fig1.pdf):
#   - data source is the new amyloid_explorer_db.tsv (829 archive
#     structures), restricted to the 775-PDB energetics cohort by
#     intersecting with Final_Structs_AE_energies.csv — vs the old metadata.txt
#     (608 / 43 / 15 / 3)
#   - column names: PDB_ID, Exp_Method etc. renamed to match the new TSV
#   - Role panel gets a Synthetic slice (~3%)
#   - Sample-origin panel gets a Recombinant bar (~62% — the visible
#     majority); the old figure dropped it
# ---------------------------------------------------------------------------

# ── Locate this script's own directory and setwd() to it ──────────────────
#
# Run-it-anywhere preamble: works whether the script is invoked via
#   `Rscript code/20260604_JS_AmylEx_Fig1.R`
# from any CWD, opened in RStudio, source()d from another R session, or
# run interactively from a plain R prompt.  Without this, `source("utils.R")`
# below silently breaks if CWD isn't `paper_figures/code/`.

.find_script_dir <- function() {
  # 1) Was the script source()d from another R session?  Walk the call
  #    stack looking for an $ofile frame variable.
  if (sys.nframe() >= 1) {
    for (i in seq_len(sys.nframe())) {
      fr <- sys.frame(i)
      ofile <- tryCatch(fr$ofile, error = function(e) NULL)
      if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
        return(dirname(normalizePath(ofile)))
      }
    }
  }
  # 2) RStudio active document?
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p)))
  }
  # 3) Rscript --file=path/to/script.R ?
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m) > 0 && file.exists(m[1])) return(dirname(normalizePath(m[1])))
  # 4) Last resort: caller's CWD.
  getwd()
}

setwd(.find_script_dir())
cat(sprintf("Working directory: %s\n", getwd()))

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggrepel)
library(stringr)
library(tibble)
library(patchwork)
library(scales)
library(readr)

# Shared palette + theme (Pathological / Functional / Synthetic; cryoEM /
# solidNMR / solutNMR; Patient / Animal / Seeded / Recombinant; etc. all
# get the SAME colour across every figure in the paper).
source("utils.R")

# Greek characters (α / β) in a PDF on macOS-from-CRAN are a rabbit hole
# (cairo needs XQuartz, ragg has no agg_pdf, etc.).  Path of least
# resistance: substitute to Latin for display only.  Underlying TSV is
# unchanged.  In Illustrator/InDesign on the final figure you can swap
# `Alpha-synuclein` → `α-synuclein` and `Abeta40/42` → `Aβ40/42` in
# under a minute (Find & Replace).

# Where the data lives + where the figure goes.  Working directory is
# already set to paper_figures/code/ by the preamble at the top of this
# script, so these relative paths are stable.
data_dir   <- file.path("..", "data")
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

# ── Load + normalise ────────────────────────────────────────────────────────

db <- readr::read_tsv(
  file.path(data_dir, "amyloid_explorer_db.tsv"),
  trim_ws = TRUE, na = c("", "NA"),
  guess_max = 1e5, show_col_types = FALSE,
  locale = locale(encoding = "UTF-8")
)

cat(sprintf("Loaded %d rows x %d cols from amyloid_explorer_db.tsv\n",
            nrow(db), ncol(db)))

# Greek-to-Latin substitution (display-only; TSV unchanged).
db <- db %>%
  mutate(
    Protein = stringr::str_replace_all(Protein,
              c("α-synuclein" = "Alpha-synuclein",
                "Aβ"          = "Abeta"))
  )

# Rename the spaced / parenthesised column names so the rest of the
# script can use bare identifiers.
db <- db %>%
  rename(
    PDB_ID        = `PDB ID`,
    Exp_Method    = `Method`,
    Role_PFS      = `Role (P/F/S)`,
    Sample_origin = `Sample origin`,
    From_Patient  = `From Patient`,
    Residues_Len  = `Residues (Length)`
  ) %>%
  mutate(
    # Coerce resolution to numeric Å
    Resolution = suppressWarnings(as.numeric(Resolution)),
    # Strip the UniProt-style organism suffix: "Homo sapiens (Human)" → "Human"
    Organism_short = case_when(
      str_detect(Organism, "Homo sapiens")                 ~ "Human",
      str_detect(Organism, "Mus musculus")                 ~ "Mouse",
      str_detect(Organism, "Rattus norvegicus")            ~ "Rat",
      str_detect(Organism, "Saccharomyces cerevisiae")     ~ "Yeast",
      str_detect(Organism, "Drosophila melanogaster")      ~ "Fly",
      str_detect(Organism, "Mesocricetus auratus")         ~ "Hamster",
      str_detect(Organism, "Gallus gallus")                ~ "Chicken",
      str_detect(Organism, "Ranoidea")                     ~ "Frog (Ranoidea)",
      str_detect(Organism, "Staphylococcus")               ~ "S. aureus",
      str_detect(Organism, "Escherichia coli")             ~ "E. coli",
      is.na(Organism) | Organism == ""                     ~ NA_character_,
      TRUE                                                  ~ Organism
    ),
    # Method normalisation: collapse RCSB long strings to publication-friendly
    # short labels.  Combined-method rows (e.g. "ELECTRON MICROSCOPY,
    # SOLID-STATE NMR") are counted under whatever appears first.
    #
    # Technique reclassification overrides:
    #   - 2BEG (Lührs et al 2005, Aβ42, PNAS) and 2NAO (Wälti/Riek 2016,
    #     Aβ42) are tagged "SOLUTION NMR" by RCSB but are in fact hybrid
    #     fibril models built from solution NMR + H/D exchange + EM-derived
    #     stack geometry + ssNMR restraints.  We classify them as solidNMR
    #     here because that's the methodology that produced the fibril
    #     structure, and because leaving them as "solutNMR" produces a
    #     misleading n=1 bar in the figure.
    #   - TODO: move this override into the Python pipeline so the public
    #     site reflects the same classification.
    Technique = case_when(
      str_detect(Exp_Method, "^ELECTRON MICROSCOPY")  ~ "cryoEM",
      str_detect(Exp_Method, "^SOLID-STATE NMR")      ~ "solidNMR",
      PDB_ID %in% c("2BEG", "2NAO")                   ~ "solidNMR",   # override (see note above)
      str_detect(Exp_Method, "^SOLUTION NMR")         ~ "solutNMR",
      str_detect(Exp_Method, "X-RAY")                 ~ "X-ray",
      TRUE                                              ~ Exp_Method
    ),
    # Role to long label
    Role_long = case_when(
      Role_PFS == "P" ~ "Pathological",
      Role_PFS == "F" ~ "Functional",
      Role_PFS == "S" ~ "Synthetic",
      TRUE             ~ "Unspecified"
    ),
    # Sample-origin to a short label for panel F
    # NOTE (2026-06-24): Nikos's official-release vocab renamed
    # recombinant_in_vitro -> in_vitro and added cell-based categories;
    # map both old and new tokens (mirrors load_energies.R).
    Origin_short = case_when(
      Sample_origin == "human_ex_vivo"                          ~ "Patient",
      Sample_origin == "animal_ex_vivo"                         ~ "Animal",
      Sample_origin == "seeded_amplified"                       ~ "Seeded",
      Sample_origin == "cellular model"                         ~ "Seeded",       # 8BGV: cellular model seeded with brain tissue
      Sample_origin == "cell extracts"                          ~ "Recombinant",  # PMEL: cell-derived functional, NOT seeded (Nikos 2026-07-14)
      Sample_origin %in% c("in_vitro", "recombinant_in_vitro")  ~ "Recombinant",
      TRUE                                                       ~ NA_character_
    )
  )

# ── Two-tier cohort accounting (Archive vs Energetics) ──────────────────────
#
# Fig1 panels show the archive composition.  Panel A info cards report
# BOTH the total archive count AND the energetics-annotated cohort count
# (the structures that survived FoldX annotation + per-residue ΔG QC).
# The donut / bar / pie panels are restricted to the energetics cohort
# so they're directly comparable with the analyses in Figs 2 onward.
#
# Single source of truth: load_energies() defines the energetics cohort
# via Final_Structs_AE_energies.csv + the per-residue ΔG sanity QC that
# excludes structures with unphysical FoldX outputs (5AEF, 2NAO under
# the current rule).  See load_energies.R for the QC details.

source("load_energies.R")
en_qc <- load_energies(keep_all_stacks = FALSE)
energy_cohort_pdbs <- unique(toupper(en_qc$PDB))

n_archive  <- nrow(db)
db_full <- db    # full archive — kept for the Panel A "Archive" card
db <- db %>% filter(toupper(PDB_ID) %in% energy_cohort_pdbs)
n_energetics <- nrow(db)
cat(sprintf(
  "Two-tier cohort: archive = %d structures, energetics = %d (dropped %d)\n",
  n_archive, n_energetics, n_archive - n_energetics
))

# Global theme — semantic palettes come from utils.R (pal_role,
# pal_technique, pal_origin, pal_species, pal_info_cards, pal_proteins).
# The Okabe-Ito raw vector is also exposed there as `okabe_ito` if you
# need a colour outside the semantic categories.
theme_set(theme_classic(base_size = 10))

# ── Tallies ─────────────────────────────────────────────────────────────────
# Two cards for the structure count (archive vs energetics-annotated)
# plus proteins.  All counts on the energetics cohort except the archive
# total.
#
# Cards DROPPED 2026-06-10:
#   - Species: 89 % of the cohort is Human; long tail dominated by
#     singletons and near-duplicate frog species.  Species composition
#     lives in the body text now.
#   - Techniques: cryoEM dominates the cohort so heavily that the count
#     is uninformative; technique breakdown is in panel G.

n_struct_archive    <- n_archive
n_struct_energetics <- n_energetics
n_protein <- n_distinct(db$Protein, na.rm = TRUE)

cat(sprintf("Totals: archive=%d, energetics=%d, %d proteins\n",
            n_struct_archive, n_struct_energetics, n_protein))

# ── Panel A — info cards ────────────────────────────────────────────────────
# Three cards: Archive structures, Energetics structures, Proteins.
# The two structure cards make the two-tier cohort accounting visible
# up front.

info_df <- tibble(
  label = factor(c("Archive\nstructures", "Energetics\nstructures",
                   "Proteins"),
                 levels = c("Archive\nstructures", "Energetics\nstructures",
                            "Proteins")),
  value = c(n_struct_archive, n_struct_energetics, n_protein),
  fill  = unname(pal_info_cards[c("Structures", "Structures", "Proteins")])
)

pA <- ggplot(info_df, aes(x = label, y = 1, fill = fill)) +
  geom_tile(width = 0.9, height = 0.85, color = "white", linewidth = 0.5, alpha = 0.9) +
  geom_text(aes(label = comma(value)), nudge_y = 0.13, fontface = "bold",
            size = 7, color = "black") +
  geom_text(aes(label = label), nudge_y = -0.20, size = 3.4, color = "black") +
  coord_cartesian(clip = "off") +
  scale_fill_identity() +
  scale_x_discrete(NULL, expand = c(0.08, 0.08)) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_void(base_size = 10) +
  theme(plot.margin = margin(t = 8, r = 8, b = 8, l = 8))

# ── Panel B — donut of top proteins ─────────────────────────────────────────

topN <- 10   # bumped from 6 — adds TMEM106B, TDP43, PrP, TAF15 etc.
             # → "Other" drops from 30% → 21% of total

prot_counts <- db %>%
  filter(!is.na(PDB_ID), !is.na(Protein), Protein != "") %>%
  count(Protein, name = "n_struct") %>%
  arrange(desc(n_struct)) %>%
  mutate(group = if_else(row_number() <= topN, Protein, "Other")) %>%
  group_by(group) %>%
  summarise(n_struct = sum(n_struct), .groups = "drop") %>%
  arrange(desc(n_struct)) %>%
  mutate(
    group = fct_reorder(group, n_struct),
    label = paste0(group, " (n=", n_struct, ")")
  ) %>%
  mutate(
    ymax = cumsum(n_struct),
    ymin = lag(ymax, default = 0),
    mid  = (ymax + ymin) / 2
  )

# Reverse for descending order in the legend
prot_counts <- prot_counts %>%
  mutate(group = factor(group, levels = rev(levels(group))))

# Use the curated 11-colour protein palette from utils.R.  "Other" always
# gets the grey at the end.
n_named <- nlevels(prot_counts$group) - 1L   # number of named proteins (everything except "Other")
fill_proteins <- setNames(
  c(pal_proteins[seq_len(n_named)], unname(okabe_ito["grey"])),
  c(setdiff(levels(prot_counts$group), "Other"), "Other")
)

pB <- ggplot(prot_counts,
             aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = group)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  scale_fill_manual(values = fill_proteins,
                    labels = paste0(prot_counts$group, " (n=", prot_counts$n_struct, ")")) +
  labs(x = NULL, y = NULL, fill = "Protein") +
  theme_void(base_size = 10) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    plot.margin       = margin(t = 6, r = 6, b = 6, l = 6)
  )

# ── Panel C — cryo-EM resolution histogram ──────────────────────────────────

cryo_df <- db %>%
  filter(Technique == "cryoEM",
         !is.na(Resolution),
         is.finite(Resolution),
         Resolution > 0, Resolution < 15)

res_med <- median(cryo_df$Resolution, na.rm = TRUE)

pC <- ggplot(cryo_df, aes(x = Resolution)) +
  # Histogram of cryo-EM resolution; bar fill = the "cryoEM" semantic colour
  # (vermillion) muted for legibility.  Stays consistent with panel G where
  # cryoEM is the dominant bar.
  geom_histogram(binwidth = 0.25,
                 fill = unname(pal_technique["cryoEM"]),
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = res_med, linetype = 2, color = "grey40") +
  annotate("text", x = res_med, y = Inf, vjust = 1.5, hjust = -0.05,
           label = sprintf("Median = %.2f Å", res_med), size = 3) +
  coord_cartesian(xlim = c(NA, 5)) +     # clip the tail at 5 Å
  labs(x = "Resolution (Å)", y = "number of structures") +
  theme(plot.title = element_blank())

# ── Panel D — Role pie (P/F/S/Unspecified) ──────────────────────────────────
# All four Role values are shown so the slices sum to the full energetics
# cohort (n = 789).  "Unspecified" (n = 76; mostly recombinant in-vitro
# α-synuclein/Tau/Aβ that the curators did not assign a P/F/S label) is a
# real curation category, NOT a dropped tail: the functional-vs-pathological
# analyses in Fig 4 / Supp Fig 4 likewise filter to Role %in% {P, F} only,
# so these structures are excluded there too.  Showing the slice here keeps
# the figure consistent with how the analysis treats them.

role_tab <- db %>%
  filter(Role_long %in% c("Pathological", "Functional", "Synthetic", "Unspecified")) %>%
  count(Role_long, name = "n") %>%
  mutate(
    Role_long = factor(Role_long,
                       levels = c("Pathological", "Functional", "Synthetic", "Unspecified")),
    frac = n / sum(n),
    lbl  = paste0(Role_long, "\nn = ", n, " (",
                  percent(frac, accuracy = 0.1), ")")
  )

# pal_role is sourced from utils.R (Pathological = orange #E69F00,
# Functional = bluish green #009E73, Synthetic = blue #0072B2).  Extend it
# locally with a neutral grey for the "Unspecified" slice.
pal_role_ext <- c(pal_role, Unspecified = "grey75")

# Pie chart of P/F/S, with labels placed by ggrepel so tiny slices
# (Synthetic at ~3%) don't get crushed.  Compute angular midpoints for
# label positioning.
role_tab <- role_tab %>%
  arrange(desc(Role_long)) %>%
  mutate(
    ymax = cumsum(frac),
    ymin = lag(ymax, default = 0),
    mid  = (ymax + ymin) / 2
  )

pD <- ggplot(role_tab,
             aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Role_long)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(-0.4, 1.6)) +
  ggrepel::geom_text_repel(
    aes(x = 1.2, y = mid, label = lbl),
    size = 2.6, segment.color = "grey40", segment.size = 0.3,
    force = 3, min.segment.length = 0,
    direction = "y", point.size = NA, max.overlaps = Inf
  ) +
  annotate("text", x = -0.4, y = 0,
           label = paste0("Total\n", sum(role_tab$n)),
           size = 2.8, fontface = "bold") +
  scale_fill_manual(values = pal_role_ext) +
  labs(x = NULL, y = NULL, fill = "Role") +
  theme_void(base_size = 10) +
  theme(legend.position = "none")

# ── Panel E — species bars (Top-2 + Other) ──────────────────────────────────

top2 <- db %>%
  filter(!is.na(Organism_short)) %>%
  count(Organism_short, name = "n") %>%
  arrange(desc(n)) %>%
  slice_head(n = 2) %>%
  pull(Organism_short)

# Keep NA hosts as an explicit "Not specified" bar so the panel sums to the
# full energetics cohort (n = 789) instead of silently dropping ~30 structures
# (mostly immunoglobulin light chains + a few recombinant/synthetic entries
# whose source organism is uncurated).  TODO: curate host for these.
species_tab <- db %>%
  mutate(group = case_when(
    is.na(Organism_short)     ~ "Not specified",
    Organism_short %in% top2  ~ Organism_short,
    TRUE                      ~ "Other"
  )) %>%
  count(group, name = "n") %>%
  mutate(
    pct = n / sum(n),
    label = paste0("n = ", comma(n), " (", percent(pct, accuracy = 0.1), ")"),
    group = factor(group, levels = c(top2, "Other", "Not specified"))
  )

pE <- ggplot(species_tab, aes(x = group, y = n, fill = group)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = label), hjust = -0.05, size = 2.9) +
  coord_flip(ylim = c(0, max(species_tab$n) * 1.2), clip = "off") +
  scale_fill_manual(values = c(pal_species, "Not specified" = "grey75"),
                    guide = "none") +    # shared palette + grey for NA host
  labs(x = NULL, y = "number of structures") +
  theme(plot.margin = margin(t = 6, r = 14, b = 6, l = 6))

# ── Panel F — sample origin bars (Patient / Animal / Seeded / Recombinant) ─

origin_tab <- db %>%
  filter(!is.na(Origin_short)) %>%
  count(Origin_short, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(
    pct = n / sum(n),
    label = paste0("n = ", comma(n), " (", percent(pct, accuracy = 0.1), ")"),
    Origin_short = factor(Origin_short,
                          levels = c("Recombinant", "Patient", "Animal", "Seeded"))
  )

# pal_origin is sourced from utils.R (Patient=yellow, Animal=sky blue,
# Seeded=green, Recombinant=reddish-purple).  Animal uses sky_blue rather
# than the deep blue #0072B2, because deep blue is reserved for "Synthetic"
# in pal_role — keeping the two distinct prevents semantic collisions
# across panels in the same figure.

pF <- ggplot(origin_tab, aes(x = Origin_short, y = n, fill = Origin_short)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = label), hjust = -0.05, size = 2.9) +
  coord_flip(ylim = c(0, max(origin_tab$n) * 1.2), clip = "off") +
  scale_fill_manual(values = pal_origin, guide = "none") +
  labs(x = NULL, y = "number of structures") +
  theme(plot.margin = margin(t = 6, r = 14, b = 6, l = 6))

# ── Panel G — technique bars ────────────────────────────────────────────────

tech_tab <- db %>%
  filter(!is.na(Technique)) %>%
  count(Technique, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(
    pct = n / sum(n),
    label = paste0("n = ", comma(n), " (", percent(pct, accuracy = 0.1), ")"),
    Technique = fct_reorder(Technique, n)
  )

# Technique fill comes from pal_technique in utils.R — so cryoEM (the
# vermillion shared with panel C's histogram), solidNMR (reddish-purple),
# solutNMR (yellow), X-ray (sky blue) all keep the same colours across
# every figure in the paper.
pG <- ggplot(tech_tab, aes(x = Technique, y = n, fill = as.character(Technique))) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = label), hjust = -0.05, size = 2.9) +
  coord_flip(ylim = c(0, max(tech_tab$n) * 1.2), clip = "off") +
  scale_fill_manual(values = pal_technique, guide = "none") +
  labs(x = NULL, y = "number of structures") +
  theme(plot.margin = margin(t = 6, r = 14, b = 6, l = 6))

# ── Assemble ────────────────────────────────────────────────────────────────
#
# Layout matches the old Fig1's three regions:
#
#   ┌──────── A (4 info cards, spans full width) ────────┐
#   │                                                     │
#   │   B (big donut)        │  C (resolution hist)       │
#   │                        │  D (Role pie)              │
#   │                        │  E (species bars)          │
#   │                                                     │
#   │   F (origin bars)      │  G (technique bars)        │
#   └─────────────────────────────────────────────────────┘
#
# patchwork design string: each letter is a 1-unit grid cell.

design <- "
AAAAAA
BBBCCC
BBBDDD
BBBEEE
FFFGGG
"

fig1 <- pA + pB + pC + pD + pE + pF + pG +
  plot_layout(design = design, heights = c(0.5, 1, 1.3, 1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

# ── Write outputs: assembled preview + individual panels ──────────────────
# See utils.R / save_panel() for the convention rationale.

# (1) Assembled patchwork — for the quick-look review only
output_pdf <- file.path(output_dir, "Fig1.pdf")
ggsave(output_pdf, fig1, width = 12, height = 10, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s (assembled preview)\n", normalizePath(output_pdf)))

# (2) Per-panel PDFs for Illustrator assembly.  Sizes are tuned to each
#     panel's natural aspect ratio (info-card strip is wide and short,
#     donut is square, bar charts are wide-and-short etc); Illustrator
#     will re-scale to final layout anyway, these are starting points.
cat("Wrote per-panel PDFs:\n")
save_panel(pA, 1, "A", width = 9.0, height = 2.5)   # info-card strip
save_panel(pB, 1, "B", width = 5.0, height = 5.0)   # donut
save_panel(pC, 1, "C", width = 5.0, height = 3.0)   # resolution histogram
save_panel(pD, 1, "D", width = 5.0, height = 4.0)   # role pie
save_panel(pE, 1, "E", width = 5.0, height = 2.5)   # species bars (3)
save_panel(pF, 1, "F", width = 5.0, height = 3.0)   # origin bars (4)
save_panel(pG, 1, "G", width = 5.0, height = 2.5)   # technique bars (3)
