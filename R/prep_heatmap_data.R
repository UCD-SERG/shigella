utils::globalVariables(c(
  ":=", ".", "sid", "unq_id", "isotype", "study_name", "serotype",
  "cohort_name", "pid", "age", "row_id", "display_label",
  "antigen", "antigen_clean", "timepoint", "timepoint_label"
))

#' Reshape + order the compiled table for the cross-reactivity heatmap
#'
#' Bangladesh (SOSAR) subset, antigens melted long, serotypes recoded (small
#' groups -> "Other"), participants ordered by serotype then age.
#'
#' @param compiled The `Compiled` sheet of the raw Shigella Excel.
#' @return A `data.table` ready for [make_heatmap()], plus a `display_label`
#'   lookup attached as attribute `"y_label_lookup"`.
#' @export
prep_heatmap_data <- function(compiled) {
  shigella <- data.table::data.table(compiled)

  data.table::setnames(
    shigella,
    old = c("sampleID", "n_ipab_MFI", "n_sf3aospbsa_MFI",
            "n_sf2aospbsa_MFI", "n_sf6ospbsa_MFI", "n_sonneiospbsa_MFI"),
    new = c("pid", "ipab", "sf3a_osp", "sf2a_osp", "sf6_osp", "sonnei_osp")
  )

  shigella[isotype == 1, sid := stringr::str_remove(unq_id, "_1:1000_IgG")]
  shigella[isotype == 2, sid := stringr::str_remove(unq_id, "_1:1000_IgA")]
  shigella[, isotype := NULL]
  data.table::setnames(shigella, old = "isotype_name", new = "isotype")

  shig <- data.table::melt(
    shigella,
    id.vars = c("pid", "sid", "study_name", "cohort", "cohort_name", "site",
                "site_name", "age", "treatment", "timepoint", "isotype"),
    measure = c("ipab", "sf3a_osp", "sf2a_osp", "sf6_osp", "sonnei_osp"),
    value.name = "result"
  )
  data.table::setnames(shig, old = "variable", new = "antigen")

  shig_bg <- shig[study_name == "SOSAR"]

  shig_bg[, serotype := data.table::fcase(
    cohort_name == "Sf2a",   "Sf2a",
    cohort_name == "Sf3a",   "Sf3a",
    cohort_name == "Sf6",    "Sf6",
    cohort_name == "sonnei", "sonnei",
    default = "Other"
  )]
  shig_bg[, serotype := factor(serotype,
                               levels = c("Sf2a", "Sf3a", "Sf6", "sonnei", "Other"))]

  participants <- shig_bg[!duplicated(pid), .(pid, serotype, age, cohort_name)]
  data.table::setorder(participants, serotype, age, pid)
  participants[, row_id := factor(pid, levels = pid)]
  participants[, display_label := ifelse(
    serotype != "Other", paste0(age, " y"),
    paste0(cohort_name, " (", age, " y)")
  )]

  shig_bg <- merge(shig_bg, participants[, .(pid, row_id, display_label)],
                   by = "pid", all.x = TRUE)

  shig_bg[, antigen_clean := data.table::fcase(
    antigen == "ipab",       "IpaB",
    antigen == "sf2a_osp",   "Sf2a",
    antigen == "sf3a_osp",   "Sf3a",
    antigen == "sf6_osp",    "Sf6",
    antigen == "sonnei_osp", "sonnei"
  )]
  shig_bg[, antigen_clean := factor(antigen_clean,
                                    levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "sonnei"))]
  shig_bg[, timepoint_label := factor(timepoint, levels = sort(unique(timepoint)))]

  attr(shig_bg, "y_label_lookup") <-
    stats::setNames(participants$display_label, as.character(participants$row_id))
  shig_bg[]
}
