#' The `codyna` Package.
#'
#' @name codyna-package
#' @description Performs analysis of complex dynamic systems with a focus on the
#' temporal unfolding of patterns, changes, and state transitions in
#' behavioral data. The package supports both time series and sequence data
#' and provides tools for the analysis and visualization of complexity,
#' pattern identification, trends, regimes, sequence typology as well as
#' early warning signals.
#'
#' @author Santtu Tikka and Mohammed Saqr
#'
"_PACKAGE"

#' Example Data on Student Engagement
#'
#' Students' engagement states (Active / Average / Disengaged)
#' throughout a whole study program. The data was generated synthetically
#' based on the article "The longitudinal association between engagement and
#' achievement varies by time, students' profiles, and achievement state:
#' A full program study". Used also in the `tna` package.
#'
#' @source \doi{10.1016/j.compedu.2023.104787}
#' @format An `stslist` object (sequence data).
"engagement"

#' Example Data on Group Regulation
#'
#' Students' regulation during collaborative learning. Students' interactions
#' were coded as:  "adapt", "cohesion", "consensus", "coregulate", "discuss",
#' "emotion", "monitor", "plan", "synthesis". Used also in the `tna` package.
#'
#' @source The data was generated synthetically.
#' @format A `data.frame` object.
"group_regulation"
