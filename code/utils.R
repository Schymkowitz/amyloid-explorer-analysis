# ---------------------------------------------------------------------------
# utils.R — shared palette + theme for the Amyloid Explorer paper figures
#
# Source this from every figure script so that semantic categories
# (Pathological / Functional / Synthetic; cryoEM / solidNMR / solutNMR;
# Patient / Animal / Seeded / Recombinant; Human / Mouse / Other) get the
# SAME colour across every figure in the paper.
#
# Usage in a figure script:
#
#     source("utils.R")     # if the script lives in code/
#     ...
#     scale_fill_manual(values = pal_role)        # P / F / S
#     scale_fill_manual(values = pal_technique)   # cryoEM / solidNMR / solutNMR
#     scale_fill_manual(values = pal_origin)      # Patient / Animal / Seeded / Recombinant
#     scale_fill_manual(values = pal_species)     # Human / Mouse / Other
#
# All hues come from the Okabe-Ito 8-colour palette — colour-blind safe,
# print-friendly, the de-facto standard for scientific figures.
# ---------------------------------------------------------------------------

# Okabe-Ito raw palette (the eight values + grey)
okabe_ito <- c(
  black          = "#000000",
  orange         = "#E69F00",
  sky_blue       = "#56B4E9",
  bluish_green   = "#009E73",
  yellow         = "#F0E442",
  blue           = "#0072B2",
  vermillion     = "#D55E00",
  reddish_purple = "#CC79A7",
  grey           = "#9E9E9E"
)

# ── Semantic palettes ──────────────────────────────────────────────────────

# Amyloid biology role — teal/orange binary scheme (2026-06-24 palette refresh).
# Functional vs Pathological is the paper's core binary comparison; teal/orange
# is colour-blind-safe, print- and projector-friendly, high contrast.
pal_role <- c(
  "Pathological" = "#D95F02",   # orange — "alert"
  "Functional"   = "#1B9E77",   # teal   — natural / growth
  "Synthetic"    = "#7570B3",   # muted purple — distinct 3rd category
  "Unspecified"  = "#9E9E9E"    # grey
)

# Experimental technique
pal_technique <- c(
  "cryoEM"   = unname(okabe_ito["vermillion"]),
  "solidNMR" = unname(okabe_ito["reddish_purple"]),
  "solutNMR" = unname(okabe_ito["yellow"]),
  "X-ray"    = unname(okabe_ito["sky_blue"]),
  "Other"    = unname(okabe_ito["grey"])
)

# Sample origin — coherent sequential blue ramp (2026-06-24 palette refresh):
# lightest = in-vitro/recombinant, darkest = direct patient ex vivo. Encodes
# "distance from the patient" as colour value; no weak yellow.
pal_origin <- c(
  "Patient"     = "#08306B",   # darkest navy — direct patient ex vivo
  "Animal"      = "#1F78B4",   # medium blue  — animal ex vivo
  "Seeded"      = "#A6CEE3",   # light blue   — seeded / amplified
  "Recombinant" = "#D9D9D9"    # light grey   — in vitro / recombinant
)

# Species (Top-2 + Other style)
pal_species <- c(
  "Human" = unname(okabe_ito["orange"]),
  "Mouse" = unname(okabe_ito["bluish_green"]),
  "Other" = unname(okabe_ito["grey"])
)

# Database-totals info cards (Fig1 panel A).  Subtle tonal variation
# rather than four hard primary colours.
pal_info_cards <- c(
  "Structures"  = "#7BAFD4",   # muted sky blue
  "Proteins"    = "#5BB28F",   # muted green
  "Species"     = "#E8C95B",   # muted yellow
  "Techniques"  = "#D89A56"    # muted orange
)

# Amino-acid biochemistry classes (Fig 3 panels C/D/E — coloring AAs
# in scatter plots by side-chain class).  Standard biochemistry-textbook
# colour conventions, chosen from Okabe-Ito where possible so the rest
# of the figure stays colour-blind safe.
pal_aa_class <- c(
  "Nonpolar"  = unname(okabe_ito["orange"]),         # hydrophobic
  "Aromatic"  = unname(okabe_ito["reddish_purple"]),
  "Polar"     = unname(okabe_ito["bluish_green"]),
  "Negative"  = unname(okabe_ito["vermillion"]),
  "Positive"  = unname(okabe_ito["blue"]),
  "Special"   = unname(okabe_ito["grey"])
)

# AA → class lookup (used by Fig 3 panels C/D/E)
aa_class <- c(
  A = "Nonpolar", V = "Nonpolar", L = "Nonpolar", I = "Nonpolar", M = "Nonpolar",
  F = "Aromatic", W = "Aromatic", Y = "Aromatic",
  P = "Special",  G = "Special",
  S = "Polar",    T = "Polar",    C = "Polar",    N = "Polar",    Q = "Polar",
  D = "Negative", E = "Negative",
  K = "Positive", R = "Positive", H = "Positive"
)

# Residue quadrant categories (Fig 3 panels H and I).  Each residue is
# scored on two axes — mean ΔG across structures (x) and SD across
# structures (y) — and assigned to one of four quadrants relative to
# the global median cut-lines:
#   core_stable            = low mean   + low SD  (stable, conserved)
#   flexible_stabilizing   = low mean   + high SD (stable, plastic)
#   rigid_destabilizing    = high mean  + low SD  (destabilizing, conserved)
#   unstable_destabilizing = high mean  + high SD (destabilizing, plastic)
#
# Pattern: cool colours for stable (low mean energy), warm colours for
# destabilizing.  Within each pair, hue distinguishes rigid vs flexible.
pal_quadrant <- c(
  "core_stable"            = unname(okabe_ito["bluish_green"]),    # stable + conserved
  "flexible_stabilizing"   = unname(okabe_ito["sky_blue"]),        # stable + plastic
  "rigid_destabilizing"    = unname(okabe_ito["reddish_purple"]),  # destabilizing + conserved
  "unstable_destabilizing" = unname(okabe_ito["vermillion"])       # destabilizing + plastic
)

# Categorical ramp for the top-N proteins donut.  Calmer than Set3's
# 12-pastel rainbow: a soft set of 11 distinguishable hues that don't
# scream.  "Other" gets the grey at the end.
pal_proteins <- c(
  "#5B9BD5",   # blue
  "#A85AB6",   # purple
  "#5BC0A8",   # teal
  "#E69F00",   # orange
  "#76A648",   # olive
  "#D38AC7",   # pink
  "#7B7DBF",   # indigo
  "#C9A86C",   # tan
  "#D55E00",   # vermillion
  "#56B4E9",   # sky blue
  "#9E9E9E"    # grey — reserved for "Other"
)

# ── Figure theme ───────────────────────────────────────────────────────────

# Consistent typography across panels.
theme_paper <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = "grey30"),
      axis.title    = ggplot2::element_text(size = base_size),
      axis.text     = ggplot2::element_text(size = base_size - 1),
      legend.title  = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text   = ggplot2::element_text(size = base_size - 2),
      legend.key.size = grid::unit(0.4, "cm"),
      plot.margin   = ggplot2::margin(t = 6, r = 8, b = 6, l = 6)
    )
}

# ── Two-track figure output ────────────────────────────────────────────────
#
# Every figure script writes TWO things:
#
#   1. an assembled patchwork preview as `output/Fig<N>.pdf`
#      → for a quick "is the data right?" look in Preview
#
#   2. one standalone PDF per panel under `output/Fig<N>_panels/`,
#      named `Fig<N>_<letter>.pdf` (e.g. Fig3_A.pdf, Fig3_B.pdf, …)
#      → these are what Joost imports into Illustrator for final
#        assembly.  Patchwork's auto-spacing is fine for previews but
#        bad for typography; per-panel PDFs let Illustrator be in
#        charge of the master layout.
#
# Use `save_panel()` for (2).  Use plain `ggsave()` for (1) (it's just
# one line in each script — not worth wrapping).

#' Save a single panel as a standalone PDF, sized for Illustrator import.
#'
#' Defaults to a 4x4-inch panel.  Override `width` / `height` per panel
#' for wide aspect ratios (heatmaps, histograms, info-card strips, etc).
#'
#' `useDingbats = FALSE` so scatter-plot point glyphs import cleanly
#' into Adobe Illustrator (Dingbats-encoded shapes sometimes corrupt).
#'
#' @param plot          ggplot or patchwork object
#' @param fig_num       integer or character — figure number (e.g. 1, 3, "Supp1")
#' @param panel_letter  character — panel id (e.g. "A", "B", "F")
#' @param width,height  numeric inches
#' @param output_dir    where to write; defaults to ../output (i.e. the
#'                      paper_figures/output/ used by every script)
#' @return invisibly, the path written.
save_panel <- function(plot, fig_num, panel_letter,
                       width  = 4,
                       height = 4,
                       output_dir = file.path("..", "output")) {
  panel_dir <- file.path(output_dir, sprintf("Fig%s_panels", fig_num))
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(panel_dir, sprintf("Fig%s_%s.pdf", fig_num, panel_letter))
  ggplot2::ggsave(path, plot,
                  width = width, height = height, units = "in",
                  device = cairo_pdf)
  message(sprintf("  panel -> %s", path))
  invisible(path)
}
