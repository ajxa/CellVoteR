#' Silently evaluate an R expression
#'
#' Evaluates an expression while suppressing all messages, warnings, and console output.
#' Useful for running noisy functions without cluttering the console.
#'
#' @param expr An R expression to evaluate.
#'
#' @return The result of the evaluated expression.
#' @keywords internal
.hush <- function(expr) {
  trash <- utils::capture.output({
    res <- suppressMessages(expr)
  })
  return(res)
}

#' Check for missing arguments in parent call
#'
#' Checks if any specific arguments are missing in the function calling this helper.
#' This looks into the `parent.frame()` to inspect the arguments of the caller.
#'
#' @param required_args A character vector of argument names to check for presence.
#'
#' @return NULL. Throws an error if any required arguments are missing.
#' @keywords internal
.check_missing_args <- function(required_args) {
  env <- parent.frame()

  is_missing <- vapply(required_args, function(arg) {
    eval(call("missing", as.name(arg)), envir = env)
  }, logical(1))

  if (any(is_missing)) {
    missing_names <- required_args[is_missing]

    cli::cli_abort(c(
      "x" = "The following required arguments are missing:",
      "*" = "{.val {missing_names}}",
      "i" = "Please provide these values in the function call."
    ))
  }
}
