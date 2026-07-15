# ---------------------------------------------------------------------------
# load_energies.R — single canonical loader for the per-residue energetics
#   table (Nikos's `Final_Structs_AE_energies.csv`) joined with the StampDB
#   public metadata (`amyloid_explorer_db.tsv`).
#
# Source this from every figure script that uses energetics so that:
#
#   - the semicolon delimiter + US decimal mark are handled in one place
#   - protein / role / sample-origin / mutation columns (which are EMPTY
#     in Nikos's CSV) get filled in from the public DB
#   - PDB ID casing is normalised (Nikos uses lower-case, our DB uses
#     mixed; the join is on upper-case)
#   - the auxiliary `Abetaa4.csv` bind_rows step from the 2025 R scripts
#     becomes unnecessary, because Abeta40/42 deposits are already in
#     the 775-structure energy cohort
#
# Usage in a figure script:
#
#     source("load_energies.R")
#     en <- load_energies()                # tidy tibble, see schema below
#     en %>% filter(Protein == "Tau") %>% ...
#
# Schema of the returned tibble (one row per (PDB, residue)):
#
#   PDB                              upper-case PDB id          (from CSV)
#   AA                               3-letter residue           (from CSV)
#   AA_1L                            1-letter residue           (from CSV)
#   nA                               residue number             (from CSV)
#   Average_Energy_Stack_1           raw avg stack energy       (from CSV)
#   Avg_SW_Energy_Stack_1            solvation-window energy    (from CSV)
#   Scaled_Average_Energy_Stack_1    z-scored vs cohort         (from CSV)
#   Scaled_Avg_SW_Energy_Stack_1     z-scored SW                (from CSV)
#   ...similarly for Stacks 2..6 (mostly NA — depend on structure)...
#   Protein                          curated short name         (from DB join)
#   Role                             P / F / S / Unspecified    (from DB join)
#   Sample_origin                    short label (Patient / …)  (from DB join)
#   Mutation                         mutation string            (from DB join)
#   From_Patient                     binary                     (from DB join)
#   Organism_short                   Human / Mouse / Other      (from DB join)
#   Technique                        cryoEM / solidNMR / …      (from DB join)
#
# All hex/colour decisions still live in utils.R; this file is data only.
# ---------------------------------------------------------------------------

# Required packages
.required <- c("dplyr", "readr", "stringr", "tibble")
for (.pkg in .required) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop(sprintf("load_energies.R: required package '%s' is not installed", .pkg))
  }
}

# ---------------------------------------------------------------------------
# Path resolution: scripts are expected to live in paper_figures/code/ and
# `source("load_energies.R")` from that directory.  Data files live in
# ../data/, populated by ../prep.sh.
# ---------------------------------------------------------------------------
.figure_data_dir <- function() {
  # Allow override (handy when sourcing from a notebook)
  if (nzchar(Sys.getenv("AE_FIGURE_DATA_DIR"))) {
    return(Sys.getenv("AE_FIGURE_DATA_DIR"))
  }
  # Otherwise, assume CWD is paper_figures/code/
  file.path("..", "data")
}

# ---------------------------------------------------------------------------
# Normalise the public DB columns so they have R-friendly bare names
# (matches the renaming the Fig1 script does — keep consistent across all
# figures).  Also derives a few short-label helper columns.
# ---------------------------------------------------------------------------
.normalise_db <- function(db) {
  dplyr::mutate(
    dplyr::rename(db,
      PDB_ID        = `PDB ID`,
      Exp_Method    = `Method`,
      Role_PFS      = `Role (P/F/S)`,
      Sample_origin = `Sample origin`,
      From_Patient  = `From Patient`,
      Residues_Len  = `Residues (Length)`
    ),
    # uppercase PDB for joining with the energy CSV
    PDB = toupper(PDB_ID),
    # Role to long label
    Role = dplyr::case_when(
      Role_PFS == "P" ~ "Pathological",
      Role_PFS == "F" ~ "Functional",
      Role_PFS == "S" ~ "Synthetic",
      TRUE             ~ "Unspecified"
    ),
    # Sample-origin short label.
    # NOTE (2026-06-23): Nikos's official-release vocab renamed
    # `recombinant_in_vitro` -> `in_vitro` (broader: includes synthetically
    # made, not only recombinant) and added cell-based categories. Both the old
    # and new tokens are mapped so the loader is robust to either DB version.
    Origin_short = dplyr::case_when(
      Sample_origin == "human_ex_vivo"                          ~ "Patient",
      Sample_origin == "animal_ex_vivo"                         ~ "Animal",
      Sample_origin == "seeded_amplified"                       ~ "Seeded",
      Sample_origin == "cellular model"                         ~ "Seeded",       # 8BGV: cellular model seeded with brain tissue
      Sample_origin == "cell extracts"                          ~ "Recombinant",  # PMEL: cell-derived functional, NOT seeded (Nikos 2026-07-14)
      Sample_origin %in% c("in_vitro", "recombinant_in_vitro")  ~ "Recombinant",
      TRUE                                                       ~ NA_character_
    ),
    # Organism normalisation
    Organism_short = dplyr::case_when(
      stringr::str_detect(Organism, "Homo sapiens")             ~ "Human",
      stringr::str_detect(Organism, "Mus musculus")             ~ "Mouse",
      stringr::str_detect(Organism, "Rattus norvegicus")        ~ "Rat",
      stringr::str_detect(Organism, "Saccharomyces cerevisiae") ~ "Yeast",
      is.na(Organism) | Organism == ""                          ~ NA_character_,
      TRUE                                                       ~ Organism
    ),
    # Method short label
    Technique = dplyr::case_when(
      stringr::str_detect(Exp_Method, "^ELECTRON MICROSCOPY") ~ "cryoEM",
      stringr::str_detect(Exp_Method, "^SOLID-STATE NMR")     ~ "solidNMR",
      stringr::str_detect(Exp_Method, "^SOLUTION NMR")        ~ "solutNMR",
      stringr::str_detect(Exp_Method, "X-RAY")                ~ "X-ray",
      TRUE                                                     ~ Exp_Method
    )
  )
}

# ---------------------------------------------------------------------------
# Read Nikos's energy CSV.  The file is semicolon-separated with US-style
# decimal `.` (NOT European `,`), so `read.csv2()` would misparse it.  We
# use readr::read_delim with explicit locale.
# ---------------------------------------------------------------------------
.read_energies <- function(path) {
  readr::read_delim(
    path,
    delim          = ";",
    locale         = readr::locale(decimal_mark = ".", grouping_mark = ""),
    na             = c("", "NA", "NaN", "nan"),
    show_col_types = FALSE,
    progress       = FALSE
  ) %>%
    dplyr::mutate(PDB = toupper(PDB)) %>%
    # Drop the empty Protein column from the CSV — we get the real one
    # from the public-DB join.
    dplyr::select(-dplyr::any_of("Protein")) %>%
    # Drop residues with unknown / non-standard AA flag
    dplyr::filter(!is.na(AA_1L), AA_1L != "X")
}

# ---------------------------------------------------------------------------
# Public entry point.
#
# Args:
#   energy_file  : filename inside data_dir (default: Final_Structs_AE_energies.csv)
#   db_file      : filename inside data_dir (default: amyloid_explorer_db.tsv)
#   data_dir     : where to find the inputs (default: ../data; respects
#                  AE_FIGURE_DATA_DIR env var)
#   keep_all_stacks: if TRUE (default), retain Stack_2..6 columns.  Set
#                  FALSE for memory-light figures that only need Stack_1.
#
# Returns: tibble keyed on (PDB, nA).  Side effect: one summary line to
# stderr so you see how many residues / PDBs got loaded.
# ---------------------------------------------------------------------------
load_energies <- function(energy_file    = "Final_Structs_AE_energies.csv",
                          db_file        = "amyloid_explorer_db.tsv",
                          data_dir       = .figure_data_dir(),
                          keep_all_stacks = TRUE,
                          apply_dG_qc    = TRUE) {

  energy_path <- file.path(data_dir, energy_file)
  db_path     <- file.path(data_dir, db_file)
  if (!file.exists(energy_path)) {
    stop(sprintf("load_energies(): cannot find %s — did you run prep.sh?", energy_path))
  }
  if (!file.exists(db_path)) {
    stop(sprintf("load_energies(): cannot find %s — did you run prep.sh?", db_path))
  }

  en <- .read_energies(energy_path)

  # ── Cohort-integrity guard (added 2026-07-13) ─────────────────────────────
  # Catches the "energy lives in a higher stack, Stack_1 empty" bug: peripheral
  # short peptides sometimes took the Stack_1 slot, pushing the main-chain
  # energies to Stack_2+, so Stack_1-only analyses silently dropped those
  # structures (8 pathological+synthetic in the 2026-07 cohort; Nikos shifted
  # them back). Warn LOUDLY and name offenders rather than let a downstream
  # filter(!is.na(Stack_1)) hide them. Never silently drop.
  {
    s1 <- "Average_Energy_Stack_1"
    hi <- grep("^Average_Energy_Stack_[2-6]$", names(en), value = TRUE)
    if (s1 %in% names(en) && length(hi)) {
      offenders <- en %>%
        dplyr::group_by(PDB) %>%
        dplyr::summarise(
          n_s1 = sum(!is.na(.data[[s1]])),
          n_hi = sum(rowSums(!is.na(dplyr::across(dplyr::all_of(hi)))) > 0),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_s1 == 0, n_hi > 0)
      if (nrow(offenders) > 0) {
        warning(sprintf(
          paste0("load_energies(): COHORT-INTEGRITY — %d structure(s) have an EMPTY ",
                 "Average_Energy_Stack_1 but energy in a higher stack; Stack_1-only ",
                 "analyses will SILENTLY DROP them: %s"),
          nrow(offenders), paste(offenders$PDB, collapse = ", ")
        ), call. = FALSE)
      } else {
        message("load_energies(): cohort-integrity OK (no empty-Stack_1 structures).")
      }
    }
  }

  db_raw <- readr::read_tsv(
    db_path,
    trim_ws        = TRUE,
    na             = c("", "NA"),
    guess_max      = 1e5,
    show_col_types = FALSE,
    locale         = readr::locale(encoding = "UTF-8")
  )

  db <- .normalise_db(db_raw)

  # Optionally drop heavy Stack_2..6 columns
  if (!keep_all_stacks) {
    drop_cols <- grep("_Stack_[2-6]$", names(en), value = TRUE)
    if (length(drop_cols)) {
      en <- dplyr::select(en, -dplyr::all_of(drop_cols))
    }
  }

  # Left-join: every energy row keeps its data, gains metadata.  Inner
  # join would silently drop residues whose PDB isn't in the public DB
  # (would happen if Nikos's cohort and our public DB are out of sync;
  # we explicitly want to surface those rows with NA biology rather than
  # hide them).
  merged <- dplyr::left_join(
    en,
    dplyr::select(db,
                  PDB,
                  Protein,
                  Role, Role_PFS,
                  Sample_origin, Origin_short,
                  From_Patient,
                  Mutation,
                  Organism, Organism_short,
                  Technique,
                  Resolution),
    by = "PDB"
  )

  # ── Per-residue ΔG sanity QC ──────────────────────────────────────────
  # FoldX outputs occasionally include physically impossible per-residue
  # energies when the input coordinates are damaged (low-resolution side
  # chains, missing atoms, bad NMR ensembles).  5AEF (5 Å Aβ42, 2015) is
  # the canonical case: I40 reports +37 kcal/mol and the structure-level
  # SD ΔG spikes to 8.3 versus a cohort typical of 1.4.  Audit of the
  # 791-PDB cohort (2026-06-10) found two structures with ≥2 residues
  # above 10 kcal/mol absolute (5AEF; 2NAO).  We flag and exclude any
  # PDB where this is true.  Threshold rationale: 99.5 % of structures
  # have max|ΔG| ≤ 10.5; values above 15 are physically unphysical for
  # a folded β-sheet residue; the AND-clause on 2-residue count is what
  # distinguishes "systematic coordinate damage" from "one anomalous
  # residue at a terminus" (which is biologically plausible).
  if (apply_dG_qc) {
    e_col <- "Average_Energy_Stack_1"
    if (e_col %in% names(merged)) {
      bad_pdbs <- merged %>%
        dplyr::filter(!is.na(.data[[e_col]])) %>%
        dplyr::group_by(PDB) %>%
        dplyr::summarise(
          max_abs = max(abs(.data[[e_col]]), na.rm = TRUE),
          n_over_10 = sum(abs(.data[[e_col]]) > 10, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::filter(max_abs > 15 | n_over_10 >= 2) %>%
        dplyr::pull(PDB)
      if (length(bad_pdbs) > 0) {
        message(sprintf(
          "load_energies(): per-residue ΔG QC excluded %d PDB(s) with unphysical energies: %s",
          length(bad_pdbs), paste(bad_pdbs, collapse = ", ")
        ))
        merged <- dplyr::filter(merged, !(PDB %in% bad_pdbs))
      }
    }
  }

  # Summary line — useful for catching e.g. "I joined the wrong DB"
  n_res   <- nrow(merged)
  n_pdb   <- dplyr::n_distinct(merged$PDB)
  n_unjoined <- dplyr::n_distinct(merged$PDB[is.na(merged$Protein)])
  message(sprintf(
    "load_energies(): %d residue rows over %d PDBs (%d PDBs in energy CSV had no metadata match)",
    n_res, n_pdb, n_unjoined
  ))

  merged
}
