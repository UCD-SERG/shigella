# ============================================================
# data.R
# Documentation for all exported data objects in the shigella package.
# ============================================================

# ── Overall population models ──────────────────────────────────────────────

#' Overall Hierarchical Model Fit: IpaB
#'
#' Fitted Bayesian hierarchical two-phase antibody kinetic model for
#' invasion plasmid antigen B (IpaB), pooling all confirmed-*Shigella*
#' cases regardless of infecting serotype.
#'
#' @format An `sr_model` tibble with columns `Iteration`, `Chain`,
#'   `Parameter`, `Iso_type`, `Stratification`, `Subject`, `value`.
#'   Key attributes:
#'   \describe{
#'     \item{`population_params`}{Tibble of `mu.par` and `prec.logy`
#'       population-level MCMC draws.}
#'     \item{`nChains`, `nIterations`}{MCMC dimensions.}
#'   }
#'
#' @source icddr,b / UCD-SERG SOSAR cohort, Dhaka 2021-2022.
#'   Mock version: 5 subjects, 20 iterations, 2 chains.
"overall_IpaB_pop_6"

#' @rdname overall_IpaB_pop_6
"overall_Sf2a_pop_6"

#' @rdname overall_IpaB_pop_6
"overall_Sf3a_pop_6"

#' @rdname overall_IpaB_pop_6
"overall_Sf6_pop_6"

#' @rdname overall_IpaB_pop_6
"overall_Sonnei_pop_6"

# ── Serotype-specific models ────────────────────────────────────────────────

#' Serotype-Specific Model Fit: S. flexneri 2a
#'
#' Hierarchical model restricted to participants infected with
#' *S. flexneri* 2a.
#' @inherit overall_IpaB_pop_6 format source
"serotype_sf2a_3"

#' @rdname serotype_sf2a_3
"serotype_sf3a_3"

#' @rdname serotype_sf2a_3
"serotype_sonnei_3"

# ── Combined S. flexneri models ────────────────────────────────────────────

#' Combined S. flexneri Model Fit: Sf2a OSP
#'
#' Pools *S. flexneri* 2a and 3a infected individuals, leveraging
#' O-antigen cross-reactivity.
#' @inherit overall_IpaB_pop_6 format source
"combined_flexneri_sf2a"

#' @rdname combined_flexneri_sf2a
"combined_flexneri_sf3a"

# ── Age-stratified IpaB models ──────────────────────────────────────────────

#' Age-Stratified IpaB Model: Children < 5 Years
#' @inherit overall_IpaB_pop_6 format source
"overall_IpaB_pop_under5"

#' Age-Stratified IpaB Model: Children >= 5 Years
#' @inherit overall_IpaB_pop_6 format source
"overall_IpaB_pop_plus5"

# ── Case data objects ───────────────────────────────────────────────────────

#' Case Data: IpaB, All Participants
#'
#' A `case_data` tibble with per-visit antibody measurements.
#'
#' @format A tibble with columns: `id`, `antigen_iso`, `timeindays`,
#'   `result`, `cohort_name`, `sid`, `isotype_name`, `visit_num`.
#'   Attributes: `id_var`, `biomarker_var`, `timeindays`, `value_var`.
#'
#' @source SOSAR cohort, Dhaka 2021-2022. Mock version: 5 subjects.
"dL_clean_Ipab_new"

#' @rdname dL_clean_Ipab_new
"dL_clean_sf2a_new"

#' @rdname dL_clean_Ipab_new
"dL_clean_sf3a_new"

#' @rdname dL_clean_Ipab_new
"dL_clean_sonnei_new"

#' @rdname dL_clean_Ipab_new
"dL_clean_sf6_new"

#' Case Data: Serotype-Restricted
#' @inherit dL_clean_Ipab_new format source
"dL_serotype_sf2a"

#' @rdname dL_serotype_sf2a
"dL_serotype_sf3a"

#' @rdname dL_serotype_sf2a
"dL_serotype_sonnei"

#' Case Data: Combined Flexneri
#' @inherit dL_clean_Ipab_new format source
"dL_combined_flexneri_sf2a_case"

#' @rdname dL_combined_flexneri_sf2a_case
"dL_combined_flexneri_sf3a_case"

#' Case Data: Age-Stratified IpaB
#' @inherit dL_clean_Ipab_new format source
"dL_clean_Ipab_new_under5"

#' @rdname dL_clean_Ipab_new_under5
"dL_clean_Ipab_new_plus5"

# ── Mock raw data (replaces Excel files) ────────────────────────────────────

#' Mock Compiled Shigella Data
#'
#' A tibble mimicking the structure of the compiled Luminex spreadsheet.
#' Used by [process_shigella_data()] and [load_all_antigens()] for
#' testing and examples.
#'
#' @format A tibble with 31 columns matching the real compiled data,
#'   including `study_name`, `sid`, `cohort_name`, `isotype_name`,
#'   `timepoint`, `Actual day`, and normalized MFI columns.
"mock_compiled_data"

#' Mock Participant Metadata
#'
#' A tibble mimicking the additional metadata spreadsheet.
#' @format A tibble with 10 columns: `sampleID`, `Gender`, `HhIncome`,
#'   `DiaBlood`, `DiaWatery`, `Fev`, `StoolNo`, `MUAC`, `IVFluidRec`,
#'   `HosDur`.
"mock_metadata"

#' Mock Duration of Diarrhea Data
#'
#' A tibble mimicking the diarrhea duration spreadsheet.
#' @format A tibble with columns `CaseID` and `DurDia_hours`.
"mock_durdia"
