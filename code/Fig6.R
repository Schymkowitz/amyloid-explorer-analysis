# ---------------------------------------------------------------------------
# 20260608_JS_AmylEx_Fig6.R — Disease-associated mutations preferentially
#   target frustrated rather than core-stabilising positions (Fig 6 in the
#   May 2026 manuscript v2.3 NL).
#
# Numbering history:
#   - revision3 (March 2026): this was Fig 7.
#   - v2.3 NL  (May  2026): renumbered to Fig 6.
#
# 3 panels:
#   A  WT vs Mutant overall stability — violin of per-structure mean ΔG,
#      grouped by WT (Mutation == NA) vs Mutant (Mutation present).  Test
#      uses MIXED-EFFECTS `mean_dG ~ is_WT + (1 | Protein)` (Type-II
#      F-test, Satterthwaite df) — same machinery as Fig 5, controlling
#      for the fact that some proteins are predominantly WT and others
#      have lots of mutants in our dataset.
#   B  Enrichment of disease mutations across four energetic contexts —
#      Core stabilising (ΔG < -0.5), Mildly stabilising (-0.5 ≤ ΔG < 0),
#      Mildly destabilising (0 ≤ ΔG < 0.5), Frustrated (ΔG ≥ 0.5).  For
#      each mutation, look up the mean WT ΔG at that residue position in
#      the corresponding protein's WT energy profile, then classify.
#      Log2 enrichment vs the background distribution of all WT residue
#      positions.  Significance by χ² goodness-of-fit against background
#      proportions.
#   C  α-synuclein WT energy profile (mean across all α-syn polymorphs)
#      with structurally resolved Parkinson's disease mutation positions
#      highlighted.  Blue ribbon = stabilising region (ΔG < 0), orange
#      ribbon = destabilising region (ΔG > 0).  Circles at mutation
#      positions, coloured by their WT energetic context.
#
# Restricted to PATHOLOGICAL structures throughout.
#
# Mutation parsing:
#   The Mutation column in amyloid_explorer_db.tsv has values like
#   "A53T" (single), "G51D, A53T" (multi), "delta19-24" or "V122del"
#   (non-point — skipped).  Position regex from the legacy script:
#     (?<=[A-Za-z])\d+(?=[A-Za-z])
#   Multi-mutation strings are split on commas before applying the regex.
#
# Source for rewrite: ANALYSIS C block in legacy code/20260324_JS_new_analyses.R
#   (lines 238-413).  Original output: fig_C_mutation_landscape.pdf in
#   old_source_material/2026/.
#
# Outputs:
#   output/Fig6.pdf                      — assembled 3-panel deliverable
#   output/Fig6_panels/Fig6_A.pdf ... C.pdf — individual panel fallbacks
# ---------------------------------------------------------------------------

# ── Run-it-anywhere preamble ──────────────────────────────────────────────
.find_script_dir <- function() {
  if (sys.nframe() >= 1) {
    for (i in seq_len(sys.nframe())) {
      fr <- sys.frame(i)
      ofile <- tryCatch(fr$ofile, error = function(e) NULL)
      if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
        return(dirname(normalizePath(ofile)))
      }
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m) > 0 && file.exists(m[1])) return(dirname(normalizePath(m[1])))
  getwd()
}
setwd(.find_script_dir())
cat(sprintf("Working directory: %s\n", getwd()))

# ── Libraries ─────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggrepel)
library(stringr)
library(tibble)
library(patchwork)
library(readr)
library(purrr)

for (pkg in c("lmerTest", "emmeans")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required: install.packages(c('lmerTest','emmeans'))"))
  }
}
library(lmerTest)
library(emmeans)

source("utils.R")
source("load_energies.R")

# ── Data ──────────────────────────────────────────────────────────────────
output_dir <- file.path("..", "output")
dir.create(output_dir, showWarnings = FALSE)

en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# Identify WT vs Mutant per PDB.  In our DB, WT structures have an NA
# (empty) Mutation field; anything non-NA is a mutant.  Restrict to
# Pathological throughout this figure.
en_patho <- en %>%
  filter(Role == "Pathological",
         !is.na(.data[[energy_col]]),
         !is.na(Protein), Protein != "") %>%
  mutate(is_WT = is.na(Mutation) | Mutation == "")

cat(sprintf("\nFig6 cohort (Pathological): %d residue rows over %d PDBs\n",
            nrow(en_patho), dplyr::n_distinct(en_patho$PDB)))

# Per-structure summary (one row per PDB) — mean ΔG + is_WT label
per_struct <- en_patho %>%
  arrange(PDB, nA) %>%
  group_by(PDB, Protein, is_WT, Mutation) %>%
  summarise(mean_dG = mean(.data[[energy_col]], na.rm = TRUE),
            .groups = "drop")

n_WT  <- sum(per_struct$is_WT)
n_Mut <- sum(!per_struct$is_WT)
cat(sprintf("  Per-structure: n_WT = %d, n_Mutant = %d\n", n_WT, n_Mut))

# ───────────────────────────────────────────────────────────────────────────
# Panel A — WT vs Mutant overall stability (mixed-effects)
# ───────────────────────────────────────────────────────────────────────────

pf_dat_A <- per_struct %>%
  mutate(WT_label = factor(if_else(is_WT, "Wild-type", "Mutant"),
                            levels = c("Wild-type", "Mutant")))

# Mixed-effects: mean_dG ~ is_WT + (1 | Protein)
m_A <- lmerTest::lmer(mean_dG ~ WT_label + (1 | Protein),
                       data = pf_dat_A, REML = TRUE)
a_A <- anova(m_A, type = "III", ddf = "Satterthwaite")
p_A <- a_A[["Pr(>F)"]][1]
cat(sprintf("\nPanel A mixed-effects: F(%g, %.1f) = %.3f, p = %.4g\n",
            a_A[["NumDF"]][1], a_A[["DenDF"]][1],
            a_A[["F value"]][1], p_A))

pal_WT_Mut <- c("Wild-type" = "#1B9E77",   # teal/orange binary scheme (2026-06-24)
                "Mutant"    = "#D95F02")

y_min_A <- min(pf_dat_A$mean_dG, na.rm = TRUE)
y_max_A <- max(pf_dat_A$mean_dG, na.rm = TRUE)
y_rng_A <- y_max_A - y_min_A
y_top_A <- y_max_A + 0.15 * y_rng_A

pA <- ggplot(pf_dat_A, aes(x = WT_label, y = mean_dG, fill = WT_label)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.16, fill = "white", color = "black",
               linewidth = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 0.7, alpha = 0.45,
              color = "grey25", stroke = 0) +
  annotate("text", x = 1.5, y = y_top_A,
           label = sprintf("Mixed-effects p = %.3g\n(WT-vs-Mut | Protein random)\nn_WT=%d  n_Mut=%d",
                            p_A, n_WT, n_Mut),
           size = 3.0, hjust = 0.5, vjust = 0, fontface = "italic") +
  scale_fill_manual(values = pal_WT_Mut, guide = "none") +
  coord_cartesian(ylim = c(y_min_A - 0.05 * y_rng_A,
                           y_top_A + 0.08 * y_rng_A),
                  clip = "off") +
  labs(x = NULL, y = expression(Mean ~ Delta * italic(G) ~ (kcal/mol))) +
  theme_paper(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 11, face = "bold"),
    plot.margin  = margin(t = 18, r = 8, b = 6, l = 6)
  )

# ───────────────────────────────────────────────────────────────────────────
# Panels B & C — shared data prep: WT energy profiles + parsed mutations
# ───────────────────────────────────────────────────────────────────────────

# Mean WT ΔG per (Protein, residue position).  Used by both Panel B (to
# classify mutation contexts) and Panel C (the α-syn profile plot).
wt_profiles <- en_patho %>%
  filter(is_WT) %>%
  group_by(Protein, nA, AA_1L) %>%
  summarise(mean_dG_WT = mean(.data[[energy_col]], na.rm = TRUE),
            n_struct   = dplyr::n_distinct(PDB),
            .groups    = "drop")

cat(sprintf("\nWT energy profiles: %d (Protein, residue) rows over %d proteins\n",
            nrow(wt_profiles), dplyr::n_distinct(wt_profiles$Protein)))

# Parse mutation positions from non-WT Mutation strings.  Handle multi-
# mutation strings like "G51D, A53T" by splitting on commas.  Drop entries
# that can't be parsed as point mutations (e.g. "delta19-24", "V122del").
mutant_info <- per_struct %>%
  filter(!is_WT, !is.na(Mutation), Mutation != "") %>%
  select(PDB, Protein, Mutation) %>%
  distinct() %>%
  mutate(mut_tokens = str_split(Mutation, ",\\s*")) %>%
  unnest(mut_tokens) %>%
  mutate(mut_token = str_trim(mut_tokens),
         mut_pos   = as.integer(str_extract(mut_token, "(?<=[A-Za-z])\\d+(?=[A-Za-z])"))) %>%
  filter(!is.na(mut_pos)) %>%
  distinct(Protein, mut_token, mut_pos)

cat(sprintf("\nParsed mutation tokens: %d unique (protein, mutation) pairs over %d proteins\n",
            nrow(mutant_info), dplyr::n_distinct(mutant_info$Protein)))

# Join mutations to the WT energy profile to get the context energy at
# each mutation site.  Drop mutations whose position isn't resolved in
# any WT structure of that protein.
mut_mapped <- mutant_info %>%
  inner_join(wt_profiles, by = c("Protein" = "Protein", "mut_pos" = "nA"))

cat(sprintf("Mutations with resolved WT context: %d (after dropping non-resolved positions)\n",
            nrow(mut_mapped)))

# ── 4-category energy context (used by B and C) ───────────────────────────
ENERGY_BREAKS <- c(-Inf, -0.5, 0, 0.5, Inf)
# Factor levels / palette keys are ASCII-only so matching is encoding-safe:
# dplyr re-encodes strings to UTF-8, which broke name-matching against
# native-encoded Unicode labels (grey bars + "NA" ticks).  The Unicode
# (ΔG ≤ ≥) lives only in the plotmath display labels below.
ENERGY_LABELS <- c("Core stabilising", "Mildly stabilising",
                   "Mildly destabilising", "Frustrated")

# Discrete RdBu, endpoints matched to the energetics gradient (#2166AC..#B2182B)
# so the 4 energetic-context bins read as the same blue->red as the ΔG heatmaps.
pal_ec <- c(
  "Core stabilising"     = "#2166AC",   # deep blue  (= energy-gradient low)
  "Mildly stabilising"   = "#92C5DE",   # light blue
  "Mildly destabilising" = "#F4A582",   # light red
  "Frustrated"           = "#B2182B"    # deep red   (= energy-gradient high)
)

# Plotmath display labels for the four contexts.  The base "pdf" device cannot
# render literal Unicode ΔG / ≤ / ≥ (they drop to "."), so the axis tick and
# legend text are rendered via plotmath (Symbol font) instead.  ec_labeller
# maps each factor level to its expression; data + palette keys stay unchanged.
ENERGY_LABELS_PM <- c(
  'atop("Core stabilising", Delta * italic(G) < -0.5)',
  'atop("Mildly stabilising", -0.5 <= Delta * italic(G) ~ "<" ~ 0)',
  'atop("Mildly destabilising", 0 <= Delta * italic(G) ~ "<" ~ 0.5)',
  'atop("Frustrated", Delta * italic(G) >= 0.5)'
)
ec_labeller <- function(x) parse(text = ENERGY_LABELS_PM[match(as.character(x), ENERGY_LABELS)])

classify_energy <- function(x) {
  cut(x, breaks = ENERGY_BREAKS, labels = ENERGY_LABELS,
      right = FALSE, include.lowest = TRUE)
}

mut_mapped <- mut_mapped %>%
  mutate(Energy_context = classify_energy(mean_dG_WT))

# ───────────────────────────────────────────────────────────────────────────
# Panel B — Enrichment of mutations across 4 energy contexts
# ───────────────────────────────────────────────────────────────────────────

# Background: distribution of all WT residue positions across all proteins.
bg_dist <- wt_profiles %>%
  mutate(Energy_context = classify_energy(mean_dG_WT)) %>%
  count(Energy_context, name = "bg_n") %>%
  mutate(bg_frac = bg_n / sum(bg_n))

# Foreground: distribution of mutation sites
mut_freq <- mut_mapped %>%
  count(Energy_context, name = "mut_n") %>%
  mutate(mut_frac = mut_n / sum(mut_n))

# Enrichment vs background.
# `classify_energy()` already returns a factor with ENERGY_LABELS as its
# levels, and `count()` on a factor preserves level order, so both
# mut_freq and bg_dist come out of the count step in the canonical
# 4-level ENERGY_LABELS order.  We rely on that order directly rather
# than re-factoring + match()ing here - the previous version of this
# block silently broke under one R/encoding combination, producing
# Energy_context = NA after the re-factor and crashing the chi-square
# call.
enrich_dat <- left_join(mut_freq, bg_dist, by = "Energy_context") %>%
  mutate(log2_enrich = log2((mut_frac + 1e-9) / (bg_frac + 1e-9))) %>%
  arrange(Energy_context)

# Defensive: assert ordering so the assumption above is checked.
stopifnot(
  "enrich_dat row count != 4 - ENERGY_LABELS mismatch" = nrow(enrich_dat) == 4,
  "bg_frac does not sum to ~1" =
    abs(sum(enrich_dat$bg_frac, na.rm = TRUE) - 1) < 1e-3
)

cat("\nMutation enrichment vs background:\n")
print(enrich_dat %>% select(Energy_context, mut_n, mut_frac, bg_n, bg_frac, log2_enrich),
      row.names = FALSE)

# χ² goodness-of-fit (4-bin omnibus, 3 df) — sledgehammer test
ct <- suppressWarnings(chisq.test(
  enrich_dat$mut_n,
  p = enrich_dat$bg_frac / sum(enrich_dat$bg_frac)   # renormalise vs rounding
))
cat(sprintf("\nχ² goodness-of-fit (4-bin omnibus): χ² = %.2f, df = %d, p = %.4f\n",
            ct$statistic, ct$parameter, ct$p.value))

# Two targeted binomial tests — these are 1-df and directly aligned with
# the figure's claim ("preferentially target frustrated rather than
# core-stabilising positions").  More powerful than the 3-df χ² because
# they don't dilute over uninteresting alternative hypotheses.
#
# Both are one-sided (alternative = "greater") because the figure is
# making a directional claim.

n_total       <- sum(enrich_dat$mut_n)

# enrich_dat is sorted in canonical ENERGY_LABELS order (1=Core stab,
# 2=Mildly stab, 3=Mildly destab, 4=Frustrated) - see the stopifnot
# guard above.  Index by position to avoid the factor-encoding trap that
# bit the chi-square block.

# (1) Frustrated (idx 4) vs the rest
n_frust       <- enrich_dat$mut_n[4]
p_frust_bg    <- enrich_dat$bg_frac[4]
bin_frust     <- binom.test(n_frust, n_total, p = p_frust_bg,
                             alternative = "greater")

# (2) Non-stabilising (idx 3,4: Mild destab + Frustrated) vs
#     Stabilising (idx 1,2: Core stab + Mild stab)
n_nonstab     <- sum(enrich_dat$mut_n[3:4])
p_nonstab_bg  <- sum(enrich_dat$bg_frac[3:4])
bin_nonstab   <- binom.test(n_nonstab, n_total, p = p_nonstab_bg,
                             alternative = "greater")

cat(sprintf("\nBinomial: Frustrated vs rest\n"))
cat(sprintf("  observed: %d/%d (%.1f%%)   |   background: %.1f%%   |   one-sided p = %.4g\n",
            n_frust, n_total, 100*n_frust/n_total,
            100*p_frust_bg, bin_frust$p.value))
cat(sprintf("\nBinomial: Non-stabilising vs Stabilising\n"))
cat(sprintf("  observed: %d/%d (%.1f%%)   |   background: %.1f%%   |   one-sided p = %.4g\n",
            n_nonstab, n_total, 100*n_nonstab/n_total,
            100*p_nonstab_bg, bin_nonstab$p.value))

y_rng_B <- max(abs(enrich_dat$log2_enrich), na.rm = TRUE)
# Three p-values stacked compactly at the top:
#   - The two 1-df binomial tests directly aligned with the figure's
#     claim (Frustrated vs rest, Non-stab vs Stab)
#   - The 3-df χ² omnibus for transparency
annot_B <- sprintf(paste(
  "Binomial (Frustrated vs rest):      p = %.3g",
  "Binomial (Non-stab vs Stab):        p = %.3g",
  "Chi-square goodness-of-fit (4-bin): p = %.3g",
  sep = "\n"),
  bin_frust$p.value, bin_nonstab$p.value, ct$p.value)

pB <- ggplot(enrich_dat, aes(x = Energy_context, y = log2_enrich,
                              fill = Energy_context)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("n=%d", mut_n),
                vjust = ifelse(log2_enrich >= 0, -0.4, 1.2)),
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = pal_ec, guide = "none") +
  scale_x_discrete(labels = ec_labeller) +
  coord_cartesian(ylim = c(-y_rng_B * 1.45, y_rng_B * 1.65), clip = "off") +
  annotate("text", x = 2.5, y = y_rng_B * 1.35,
           label = annot_B,
           size = 2.8, hjust = 0.5, vjust = 0,
           fontface = "italic", lineheight = 1.0, family = "mono") +
  labs(x = "Energetic context of mutation site (WT profile)",
       y = expression(log[2]~enrichment~vs~background)) +
  theme_paper(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8.5),
    plot.margin = margin(t = 22, r = 8, b = 6, l = 6)
  )

# ───────────────────────────────────────────────────────────────────────────
# Panel C — α-synuclein WT profile with mutation positions highlighted
# ───────────────────────────────────────────────────────────────────────────
#
# Lead protein for this analysis is α-synuclein (most studied disease
# mutations + plenty of WT structures in the cohort).  Same recipe will
# work for Tau, TTR, IAPP etc — change ASYN_NAME below to swap.

ASYN_NAME <- "α-synuclein"

# IMPORTANT: filtering on Protein == "α-synuclein" with `==` silently
# returns ZERO rows under one R/encoding combination because the
# in-source Greek α has Encoding("unknown") while the DB-loaded string
# is marked Encoding("UTF-8") - byte-identical but `==` returns FALSE.
# (Same trap that bit SuppFig3.)  Use grepl on the unique substring
# "synuclein" instead, and assert non-empty before continuing so any
# future regression is loud.
is_asyn <- function(p) grepl("synuclein", p, ignore.case = TRUE)

asyn_wt   <- wt_profiles %>% filter(is_asyn(Protein)) %>% arrange(nA)
asyn_muts <- mut_mapped  %>% filter(is_asyn(Protein))

cat(sprintf("\nPanel C (%s): %d resolved WT residues, %d resolved mutations\n",
            ASYN_NAME, nrow(asyn_wt), nrow(asyn_muts)))

stopifnot(
  "Panel C: no α-synuclein WT residues resolved - check Protein labelling / Encoding" =
    nrow(asyn_wt) > 0
)

if (nrow(asyn_wt) == 0) {
  warning(sprintf("No WT structures resolved for %s — skipping Panel C", ASYN_NAME))
  pC <- ggplot() + theme_void() +
    labs(title = sprintf("Panel C skipped: no WT %s residues resolved", ASYN_NAME))
} else {
  # Drop unused factor levels so the legend only shows contexts that
  # actually contain a mutation
  asyn_muts <- asyn_muts %>%
    mutate(Energy_context = droplevels(Energy_context))

  pC <- ggplot(asyn_wt, aes(x = nA, y = mean_dG_WT)) +
    # Blue ribbon for stabilising regions (ΔG < 0)
    geom_ribbon(aes(ymin = pmin(mean_dG_WT, 0), ymax = 0),
                fill = unname(okabe_ito["sky_blue"]), alpha = 0.45) +
    # Orange ribbon for destabilising regions (ΔG > 0)
    geom_ribbon(aes(ymin = 0, ymax = pmax(mean_dG_WT, 0)),
                fill = unname(okabe_ito["orange"]), alpha = 0.45) +
    geom_line(linewidth = 0.6, color = "black") +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed",
               color = "grey40") +
    geom_point(data = asyn_muts,
               aes(x = mut_pos, y = mean_dG_WT, fill = Energy_context),
               shape = 21, size = 4, color = "black", stroke = 0.8) +
    ggrepel::geom_text_repel(data = asyn_muts,
                              aes(x = mut_pos, y = mean_dG_WT,
                                  label = mut_token),
                              size = 3.2, segment.size = 0.3,
                              nudge_y = 0.5, force = 3,
                              min.segment.length = 0,
                              max.overlaps = Inf,
                              fontface = "bold") +
    scale_fill_manual(values = pal_ec, name = "Energetic context",
                      labels = ec_labeller, drop = FALSE) +
    labs(x = "Residue position",
         y = expression(Mean ~ WT ~ Delta * italic(G) ~ (kcal/mol))) +
    theme_paper(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 8),
          plot.margin     = margin(t = 6, r = 8, b = 6, l = 6))
}

# ── Assemble + save ───────────────────────────────────────────────────────
fig6 <- ((pA | pB) / pC) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

output_pdf <- file.path(output_dir, "Fig6.pdf")
ggsave(output_pdf, fig6, width = 14, height = 11, units = "in",
       device = cairo_pdf)
cat(sprintf("\nWrote %s\n", normalizePath(output_pdf)))

cat("Wrote per-panel PDFs:\n")
save_panel(pA, 6, "A", width = 4.5, height = 5.0)
save_panel(pB, 6, "B", width = 6.5, height = 5.0)
save_panel(pC, 6, "C", width = 13.0, height = 6.0)
