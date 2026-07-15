# ---------------------------------------------------------------------------
# 20260610_JS_AmylEx_robustness_audit.R - Bootstrap robustness audit
#   for every claim currently slated for the manuscript.  Random 80 %
#   protein subsamples, 1000 iterations (200-500 for the slow mixed-
#   effects tests).  Reports effect-size CI, % iterations with the same
#   direction as the full-cohort estimate, % iterations clearing α = 0.05,
#   and a four-bin classification {robust / stable / borderline / fragile}.
#
# Set rubric BEFORE looking at results (per Joost's principle from the
# discussion):
#   robust     : direction-match > 95 % AND effect-size 95 % CI excludes 0
#   stable     : direction-match > 90 %, CI brushes 0 only narrowly
#   borderline : direction-match 50-90 %; p flickers across α
#   fragile    : direction-match < 50 % OR direction flips
#
# Output: data/robustness_audit.tsv (one row per claim).
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

# ── Libraries ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(lmerTest)
})

source("utils.R")
source("load_energies.R")

cat(sprintf("Working directory: %s\n", getwd()))

# ═════════════════════════════════════════════════════════════════════════
# DATA LOAD (once, captured by closure into each test fn)
# ═════════════════════════════════════════════════════════════════════════
en <- load_energies(keep_all_stacks = FALSE)
energy_col <- "Average_Energy_Stack_1"

# Per-PDB metrics for Fig 4 architecture + Fig 5 mixed-effects
per_struct <- en %>%
  filter(!is.na(.data[[energy_col]])) %>%
  arrange(PDB, nA) %>%
  group_by(PDB, Role, Protein, Origin_short, Resolution) %>%
  summarise(
    mean_dG = mean(.data[[energy_col]], na.rm = TRUE),
    frac_stab = mean(.data[[energy_col]] < 0, na.rm = TRUE),
    sd_dG = stats::sd(.data[[energy_col]], na.rm = TRUE),
    mean_run_stab = {
      vals <- .data[[energy_col]]
      r <- rle(vals < 0)
      sr <- r$lengths[r$values == TRUE]
      if (length(sr) > 0) mean(sr) else NA_real_
    },
    .groups = "drop"
  )

# Surface composition
surf_raw <- read_tsv(
  file.path("..", "data", "amyloid_surface_composition_per_residue.tsv"),
  show_col_types = FALSE, progress = FALSE
) %>%
  mutate(PDB = sub("^AMEX_", "", pdb)) %>%
  inner_join(en %>% distinct(PDB, Protein, Role), by = "PDB") %>%
  filter(exposed == 1, aa1 %in% LETTERS)

# AAindex scales
props_raw <- read_tsv(file.path("..", "data", "physprops.txt"),
                     show_col_types = FALSE, progress = FALSE)
names(props_raw)[1] <- "Scale"

# All proteins
all_proteins <- unique(en$Protein[!is.na(en$Protein) & en$Protein != ""])
n_prot <- length(all_proteins)
cat(sprintf("Cohort: %d unique proteins, %d PDBs\n",
            n_prot, n_distinct(en$PDB)))

# Mutation parsing (for Fig 6)
mut_tokens <- en %>%
  filter(Role == "Pathological",
         !is.na(Mutation), Mutation != "",
         !(toupper(Mutation) %in% c("WT", "NONE", "NA"))) %>%
  distinct(PDB, Protein, Mutation) %>%
  rowwise() %>%
  mutate(tokens = list(str_extract_all(Mutation, "[A-Z][0-9]+[A-Z]")[[1]])) %>%
  ungroup() %>%
  unnest(tokens) %>%
  mutate(pos = as.integer(str_extract(tokens, "[0-9]+"))) %>%
  filter(!is.na(pos)) %>%
  distinct(Protein, pos, tokens)

# ═════════════════════════════════════════════════════════════════════════
# BOOTSTRAP UTILITY
# ═════════════════════════════════════════════════════════════════════════
bootstrap_finding <- function(test_fn, n_iter = 1000,
                              fraction = 0.8, seed = 42) {
  set.seed(seed)
  n_keep <- floor(fraction * n_prot)
  iters <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    kept <- sample(all_proteins, n_keep)
    res <- tryCatch(
      test_fn(kept),
      error = function(e) tibble(effect = NA_real_, p = NA_real_,
                                  err = conditionMessage(e))
    )
    res$iter <- i
    iters[[i]] <- res
  }
  bind_rows(iters)
}

summarise_robustness <- function(boot_df, full_effect, full_p, label,
                                  alpha = 0.05) {
  ok <- !is.na(boot_df$effect) & !is.na(boot_df$p)
  d <- boot_df[ok, ]
  if (nrow(d) == 0) {
    return(tibble(label = label, n_ok = 0, n_failed = nrow(boot_df),
                  err = "all iterations failed"))
  }
  dir_match <- mean(sign(d$effect) == sign(full_effect))
  ci_lo <- as.numeric(quantile(d$effect, 0.025))
  ci_hi <- as.numeric(quantile(d$effect, 0.975))
  ci_excludes_zero <- (ci_lo > 0 && ci_hi > 0) || (ci_lo < 0 && ci_hi < 0)
  pct_sig <- mean(d$p < alpha)

  class <- case_when(
    dir_match > 0.95 & ci_excludes_zero ~ "robust",
    dir_match > 0.90 ~ "stable",
    dir_match > 0.50 ~ "borderline",
    TRUE             ~ "fragile"
  )

  tibble(
    label         = label,
    full_effect   = full_effect,
    full_p        = full_p,
    n_ok          = nrow(d),
    n_failed      = nrow(boot_df) - nrow(d),
    median_effect = median(d$effect),
    ci025         = ci_lo,
    ci975         = ci_hi,
    pct_dir_match = dir_match,
    pct_sig       = pct_sig,
    class         = class
  )
}

# Helper: extract per-AA values for an AAindex scale.
get_scale_values <- function(scale_name) {
  row <- props_raw %>% filter(Scale == scale_name)
  if (nrow(row) == 0) return(NULL)
  v <- as.numeric(unlist(row[1, -1]))
  names(v) <- names(props_raw)[-1]
  v
}

# ═════════════════════════════════════════════════════════════════════════
# TEST FUNCTIONS
# ═════════════════════════════════════════════════════════════════════════

# --- FIG 3 scale correlations: per-AA mean ΔG vs scale values ------------
make_scale_corr_test <- function(scale_name) {
  s_vec <- get_scale_values(scale_name)
  function(prot_subset) {
    aa_energy <- en %>%
      filter(Protein %in% prot_subset) %>%
      mutate(AA_1L = toupper(AA_1L)) %>%
      filter(AA_1L %in% LETTERS, !is.na(.data[[energy_col]])) %>%
      group_by(AA_1L) %>%
      summarise(mean_energy = mean(.data[[energy_col]], na.rm = TRUE),
                .groups = "drop")
    overlap <- intersect(aa_energy$AA_1L, names(s_vec))
    if (length(overlap) < 5) return(tibble(effect = NA, p = NA))
    e_vec <- aa_energy$mean_energy[match(overlap, aa_energy$AA_1L)]
    sv <- s_vec[overlap]
    if (sd(sv, na.rm = TRUE) == 0) return(tibble(effect = NA, p = NA))
    ct <- suppressWarnings(cor.test(e_vec, sv))
    tibble(effect = unname(ct$estimate), p = ct$p.value, n_aa = length(overlap))
  }
}

# --- FIG 4 architecture (F vs P): MW U per metric, effect = Cohen's d ----
cohen_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  pooled <- sqrt(((nx - 1) * var(x, na.rm = TRUE) +
                  (ny - 1) * var(y, na.rm = TRUE)) / max(nx + ny - 2, 1))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled
}

make_arch_test <- function(metric) {
  function(prot_subset) {
    d <- per_struct %>%
      filter(Protein %in% prot_subset,
             Role %in% c("Functional", "Pathological"),
             !is.na(.data[[metric]]))
    f <- d[[metric]][d$Role == "Functional"]
    p <- d[[metric]][d$Role == "Pathological"]
    if (length(f) < 3 || length(p) < 3)
      return(tibble(effect = NA, p = NA))
    wt <- suppressWarnings(wilcox.test(f, p, exact = FALSE))
    tibble(effect = cohen_d(f, p), p = wt$p.value,
           n_F = length(f), n_P = length(p))
  }
}

# --- FIG 4 surface composition (F vs P) on a specific scale --------------
make_surf_scale_test <- function(scale_name) {
  s_vec <- get_scale_values(scale_name)
  function(prot_subset) {
    sub <- surf_raw %>%
      filter(Protein %in% prot_subset,
             Role %in% c("Functional", "Pathological"),
             aa1 %in% names(s_vec))
    if (nrow(sub) == 0) return(tibble(effect = NA, p = NA))
    sub$sv <- s_vec[sub$aa1]
    sub <- sub %>% filter(!is.na(sv))
    pdb_means <- sub %>%
      group_by(PDB, Protein, Role) %>%
      summarise(surf_mean = mean(sv, na.rm = TRUE), .groups = "drop")
    prot_means <- pdb_means %>%
      group_by(Protein, Role) %>%
      summarise(pm = mean(surf_mean, na.rm = TRUE), .groups = "drop")
    f <- prot_means$pm[prot_means$Role == "Functional"]
    p <- prot_means$pm[prot_means$Role == "Pathological"]
    if (length(f) < 3 || length(p) < 3)
      return(tibble(effect = NA, p = NA))
    wt <- suppressWarnings(wilcox.test(f, p, exact = FALSE))
    tibble(effect = cohen_d(f, p), p = wt$p.value,
           n_F_prot = length(f), n_P_prot = length(p))
  }
}

# --- FIG 5 mixed-effects Origin omnibus (Resolution-controlled) -----------
make_fig5_omnibus_test <- function(metric) {
  function(prot_subset) {
    d <- per_struct %>%
      filter(Protein %in% prot_subset,
             Role == "Pathological",
             Origin_short %in% c("Recombinant", "Seeded", "Animal", "Patient"),
             !is.na(.data[[metric]]), !is.na(Resolution))
    if (length(unique(d$Protein)) < 8) return(tibble(effect = NA, p = NA))
    if (length(unique(d$Origin_short)) < 4) return(tibble(effect = NA, p = NA))
    d$Origin_short <- factor(d$Origin_short,
                              levels = c("Recombinant", "Seeded",
                                         "Animal", "Patient"))
    fit <- tryCatch(
      lmerTest::lmer(
        as.formula(sprintf("%s ~ Origin_short + Resolution + (1|Protein)",
                           metric)),
        data = d, REML = TRUE
      ),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (is.null(fit)) return(tibble(effect = NA, p = NA))
    aov <- tryCatch(anova(fit, type = "III"),
                     error = function(e) NULL)
    if (is.null(aov) || !("Origin_short" %in% rownames(aov)))
      return(tibble(effect = NA, p = NA))
    p_origin <- aov["Origin_short", "Pr(>F)"]
    f_origin <- aov["Origin_short", "F value"]
    tibble(effect = f_origin, p = p_origin, n_prot = length(unique(d$Protein)))
  }
}

# --- FIG 6 cohort-wide binomial (non-stab vs stab) and χ² ----------------
make_fig6_binomial_test <- function() {
  function(prot_subset) {
    # WT profiles for kept proteins
    wt <- en %>%
      filter(Protein %in% prot_subset,
             (is.na(Mutation) | Mutation == "" |
                toupper(Mutation) %in% c("WT", "NONE", "NA"))) %>%
      group_by(Protein, nA) %>%
      summarise(mean_dG_WT = mean(.data[[energy_col]], na.rm = TRUE),
                .groups = "drop") %>%
      mutate(ctx = cut(mean_dG_WT,
                       breaks = c(-Inf, -0.5, 0, 0.5, Inf),
                       labels = 1:4, right = FALSE,
                       include.lowest = TRUE))
    mut_in_subset <- mut_tokens %>%
      filter(Protein %in% prot_subset) %>%
      inner_join(wt, by = c("Protein" = "Protein", "pos" = "nA")) %>%
      filter(!is.na(ctx))
    if (nrow(mut_in_subset) < 10) return(tibble(effect = NA, p = NA))
    mut_n <- as.numeric(table(factor(mut_in_subset$ctx, levels = 1:4)))
    bg_n  <- as.numeric(table(factor(wt$ctx, levels = 1:4)))
    if (any(bg_n == 0)) return(tibble(effect = NA, p = NA))
    bg_frac <- bg_n / sum(bg_n)
    n_total <- sum(mut_n)
    n_nonstab <- sum(mut_n[3:4])
    p_nonstab_bg <- sum(bg_frac[3:4])
    bin <- suppressWarnings(binom.test(n_nonstab, n_total,
                                       p = p_nonstab_bg,
                                       alternative = "greater"))
    enrich <- log2((n_nonstab / n_total + 1e-9) / (p_nonstab_bg + 1e-9))
    tibble(effect = enrich, p = bin$p.value, n_mut = n_total)
  }
}

# --- SUPP FIG 3: per-AA median ΔG correlation full vs tau+α-syn-excluded -
make_suppfig3_test <- function() {
  is_dom <- function(p) grepl("^Tau$", p) |
                         grepl("synuclein", p, ignore.case = TRUE)
  function(prot_subset) {
    sub <- en %>%
      filter(Protein %in% prot_subset) %>%
      mutate(AA_1L = toupper(AA_1L)) %>%
      filter(AA_1L %in% LETTERS, !is.na(.data[[energy_col]]))
    full_med <- sub %>%
      group_by(AA_1L) %>%
      summarise(m = median(.data[[energy_col]], na.rm = TRUE),
                .groups = "drop")
    red_med <- sub %>% filter(!is_dom(Protein)) %>%
      group_by(AA_1L) %>%
      summarise(m = median(.data[[energy_col]], na.rm = TRUE),
                .groups = "drop")
    j <- inner_join(full_med, red_med, by = "AA_1L",
                    suffix = c("_full", "_red"))
    if (nrow(j) < 5) return(tibble(effect = NA, p = NA))
    ct <- suppressWarnings(cor.test(j$m_full, j$m_red))
    tibble(effect = unname(ct$estimate), p = ct$p.value, n_aa = nrow(j))
  }
}

# ═════════════════════════════════════════════════════════════════════════
# REGISTRY OF CLAIMS TO AUDIT
# ═════════════════════════════════════════════════════════════════════════
# Per-claim: a (label, test_fn, n_iter, full_effect, full_p) tuple.
# Run once on the full cohort first to get the headline numbers, then
# bootstrap.

run_claim <- function(label, test_fn, n_iter = 1000) {
  cat(sprintf("\n[%s] full-cohort first...\n", label))
  full <- test_fn(all_proteins)
  cat(sprintf("  full effect = %s, p = %s\n",
              format(full$effect %||% NA, digits = 3),
              format(full$p %||% NA, digits = 3)))
  cat(sprintf("  bootstrapping (%d iter)...\n", n_iter))
  t0 <- Sys.time()
  boot <- bootstrap_finding(test_fn, n_iter = n_iter)
  cat(sprintf("  done in %.1f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  summarise_robustness(boot, full$effect, full$p, label)
}

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# ═════════════════════════════════════════════════════════════════════════
# RUN ALL CLAIMS
# ═════════════════════════════════════════════════════════════════════════
results <- list(
  # Fig 3 — residue grammar
  run_claim("Fig3_GOLD730101",       make_scale_corr_test("GOLD730101"), 1000),
  run_claim("Fig3_GRAR740102",       make_scale_corr_test("GRAR740102"), 1000),
  run_claim("Fig3_Aggrescan",        make_scale_corr_test("Aggrescan"),  1000),

  # Fig 4 architecture (F vs P)
  run_claim("Fig4_arch_meanDG",      make_arch_test("mean_dG"),       1000),
  run_claim("Fig4_arch_fracStab",    make_arch_test("frac_stab"),     1000),
  run_claim("Fig4_arch_meanStabRun", make_arch_test("mean_run_stab"), 1000),
  run_claim("Fig4_arch_sdDG",        make_arch_test("sd_dG"),         1000),

  # Fig 4 surface composition
  run_claim("Fig4_surf_GOLD730101",  make_surf_scale_test("GOLD730101"), 1000),
  run_claim("Fig4_surf_HOPT810101",  make_surf_scale_test("HOPT810101"), 1000),

  # Fig 5 mixed-effects Origin
  run_claim("Fig5_fracStab_Origin",    make_fig5_omnibus_test("frac_stab"),    300),
  run_claim("Fig5_sdDG_Origin",        make_fig5_omnibus_test("sd_dG"),        300),
  run_claim("Fig5_meanDG_Origin",      make_fig5_omnibus_test("mean_dG"),      300),
  run_claim("Fig5_meanStabRun_Origin", make_fig5_omnibus_test("mean_run_stab"),300),

  # Fig 6 — mutation enrichment
  run_claim("Fig6_binomial_nonstab", make_fig6_binomial_test(), 1000),

  # Supp Fig 3 — tau/α-syn-excluded robustness
  run_claim("SuppFig3_perAA_corr",   make_suppfig3_test(),      1000)
)

results_tbl <- bind_rows(results)

cat("\n\n=== ROBUSTNESS AUDIT SUMMARY ===\n")
results_tbl %>%
  select(label, full_effect, full_p, median_effect,
         ci025, ci975, pct_dir_match, pct_sig, class) %>%
  print(n = Inf, width = Inf)

out <- file.path("..", "data", "robustness_audit.tsv")
write_tsv(results_tbl, out)
cat(sprintf("\nWrote %s\n", normalizePath(out)))
