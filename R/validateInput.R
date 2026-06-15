# Shared, actionable input validation for all DNB methods.
# One place defines what a valid call looks like, so every method fails the
# same way -- a clear message telling the user what is wrong and how to fix it,
# never a cryptic error bubbling up from cor()/hclust()/subscripting.

#' @title Validate the two-column state table
#' @param state The user-supplied state data.frame.
#' @return Invisibly \code{TRUE}; stops with an actionable error otherwise.
#' @noRd
checkStateTwoCol <- function(state) {
    if (ncol(state) != 2) {
        stop("'state' must have exactly two columns: sample IDs in the first, ",
             "group/time labels in the second.", call. = FALSE)
    }
    invisible(TRUE)
}

#' @title Ensure every state level has at least one sample
#' @param ddl A list of per-state expression matrices (state -> genes x samples).
#' @return Invisibly \code{TRUE}; stops listing the empty states otherwise.
#' @noRd
checkNoEmptyState <- function(ddl) {
    empty <- names(ddl)[vapply(ddl, NCOL, integer(1)) == 0]
    if (length(empty) > 0) {
        stop("no samples found for state(s): ", paste(empty, collapse = ", "),
             ".\n  Check that 'state.levels' matches the values in the state ",
             "column,\n  and that sample IDs match colnames(expr).",
             call. = FALSE)
    }
    invisible(TRUE)
}

#' @title Ensure a reference state named "ref" is present
#' @param state.levels The user-supplied state-level vector.
#' @return Invisibly \code{TRUE}; stops with an actionable error otherwise.
#' @noRd
checkRefState <- function(state.levels) {
    if (!"ref" %in% state.levels) {
        stop("'state.levels' must include \"ref\" as the reference state.\n",
             "  Got: ", paste(state.levels, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

#' @title Ensure a count is non-zero after filtering
#' @param n The count to check (e.g. number of edges or case samples).
#' @param what A short noun describing what was counted, used in the message.
#' @param hint Optional remediation hint appended to the message.
#' @return Invisibly \code{TRUE}; stops with an actionable error otherwise.
#' @noRd
checkNonEmpty <- function(n, what, hint = "") {
    if (n == 0) {
        stop("no ", what, " left after filtering.",
             if (nzchar(hint)) paste0("\n  ", hint) else "", call. = FALSE)
    }
    invisible(TRUE)
}
