#' Print a formatted box of parameters
#'
#' prints a stylized CLI box containing a list of parameters and their values
#' to the console. It is designed for logging configuration settings at the
#' start of a script or process.
#'
#' @param param_list A named list. The names of the list elements are used as
#'   labels, and the values are displayed in blue.
#'
#' @return No return value, called for side effects (printing to console).
#'
#'
#' @examples
#' \dontrun{
#'   config <- list(
#'     date = "2023-10-27",
#'     model = "XGBoost",
#'     iterations = 100,
#'     verbose = TRUE
#'   )
#'   params_box(config)
#' }
#'
#' @export
params_box <- function(param_list) {

  max_w <- max(nchar(names(param_list))) + 10

  txt <- cli::ansi_columns(
    text = paste0(names(param_list), ": ", cli::col_blue(param_list)),
    width = max_w,
    fill = "rows",
    max_cols = 1,
    align = "left"
  )

  print(
    cli::boxx(
      label = txt, padding = c(0, 2, 1, 0),
      header = cli::col_br_white("Parameters")
    )
  )

}




#' Print a styled H1 header with optional contextual information
#'
#' This function prints a double-line horizontal rule to the console, acting as a
#' major section header. It supports a main text label and an optional, italicized
#' "info" label (e.g., for Sample IDs or timestamps), allowing for distinct
#' coloring of the line, the main text, and the info text.
#'
#' @param text Character string. The primary header text.
#' @param info Character string (optional). Secondary text to display in
#'   parentheses next to the main text. This text is automatically italicized.
#'   Defaults to \code{NULL}.
#' @param head_color Character string. Color of the horizontal rule lines.
#'   Defaults to \code{"cyan"}.
#' @param text_color Character string. Color of the main text.
#'   Defaults to \code{"white"}.
#' @param info_color Character string. Color of the secondary info text.
#'   Defaults to \code{"grey60"}.
#'
#' @return No return value; prints a formatted message to the console.
#'
#' @details
#' The function combines `cli::cli_rule` for the structure with inline ANSI
#' styling for the labels. The `info` argument, if provided, is formatted as
#' `(info)`, colored according to `info_color`, and italicized using
#' `cli::style_italic`.
#'
#' @importFrom cli cli_div cli_rule cli_end make_ansi_style style_italic
#'
#' @examples
#' \dontrun{
#'   # Standard usage
#'   print_h1("Data Ingestion")
#'
#'   # With context info (e.g. Sample ID)
#'   # Output: == CellVoterR (Sample_001) ==
#'   print_h1("CellVoterR", info = "Sample_001")
#'
#'   # Custom colors
#'   print_h1("Analysis", info = "Batch A", head_color = "red", info_color = "yellow")
#' }
#' @export
print_h1 <- function(text,
                     info = NULL,
                     head_color = "cyan",
                     text_color = "cyan",
                     info_color = "grey60") {


  main_style <- cli::make_ansi_style(text_color)
  info_style <- cli::make_ansi_style(info_color)

  if (!is.null(info) && info != "") {
    info_text <- paste0("(", info, ")")
    formatted_info <- cli::style_italic(info_style(info_text))
    label <- paste0(main_style(text), " ", formatted_info)
  } else {
    label <- main_style(text)
  }

  div_id <- cli::cli_div(
    theme = list(
      rule = list(
        color = head_color,
        "font-weight" = "bold",
        "line-type" = "double"
      )
    )
  )

  cli::cli_rule(label)
  cli::cli_end(div_id)
}



#' Print a styled H2 header
#'
#' Prints a single-line horizontal rule with a split layout. The main text appears
#' on the left (bold, colored), and optional info appears on the right (italic, grey).
#'
#' @param text Character string. The text for the left side of the header.
#' @param info Character string (Optional). Text for the right side, displayed in brackets.
#' @param line_color Character string. Color of the horizontal line. Defaults to "cyan".
#' @param text_color Character string. Color of the left-side text. Defaults to "cyan".
#' @param info_color Character string. Color of the right-side info text. Defaults to "grey60".
#'
#' @export
print_h2 <- function(text,
                     info = NULL,
                     line_color = "cyan",
                     text_color = "cyan",
                     info_color = "grey80") {

  left_style  <- cli::make_ansi_style(text_color)
  right_style <- cli::make_ansi_style(info_color)

  left_label <- cli::style_bold(left_style(text))

  right_label <- ""
  if (!is.null(info) && info != "") {
    raw_info <- paste0("(", info, ")")
    right_label <- cli::style_italic(right_style(raw_info))
  }

  div_id <- cli::cli_div(
    theme = list(
      rule = list(
        color = line_color,
        "line-type" = "single"
      )
    )
  )

  escape_cli_braces <- function(x) gsub("\\{", "{{", gsub("\\}", "}}", x))

  cli::cli_rule(
    left = escape_cli_braces(left_label),
    right = escape_cli_braces(right_label)
    )

  cli::cli_end(div_id)
}





#' Print a status bar with icons
#'
#' Prints a summary line with two statuses separated by a vertical bar.
#' Success items get a Green Tick; Failure/Unknown items get a Red Cross.
#'
#' @param n_pass Number (or percentage) for the "good" category (e.g., labelled cells).
#' @param n_fail Number (or percentage) for the "bad" category (e.g., unknown cells).
#' @param label_text Character. Label for the good category. Defaults to "labelled".
#' @param unknown_text Character. Label for the bad category. Defaults to "unknown".
#' @param use_props Logical. If TRUE, appends a "%" sign to the numbers. Defaults to FALSE.
#'
#' @export
print_status_bar <- function(n_pass, n_fail,
                             label_text = "Labelled",
                             unknown_text = "Unknown",
                             use_props = FALSE) {


  color_pass <- cli::make_ansi_style("green")
  color_fail <- cli::make_ansi_style("red")
  color_pipe <- cli::make_ansi_style("white")

  val_pass <- if (use_props) paste0(n_pass, "%") else format(n_pass, big.mark = ",")
  val_fail <- if (use_props) paste0(n_fail, "%") else format(n_fail, big.mark = ",")

  str_pass <- paste0(val_pass, " ", label_text, " ", cli::symbol$tick)
  styled_pass <- color_pass(str_pass)

  str_fail <- paste0(val_fail, " ", unknown_text, " ", cli::symbol$cross)
  styled_fail <- color_fail(str_fail)

  sep <- color_pipe(" | ")

  cli::cli_text("")
  cli::cli_text(paste0(styled_pass, sep, styled_fail))
}




#' Print an inline vector summary
#'
#' Prints a list of labeled counts separated by vertical pipes.
#' Format: "label 1 (count 1) | label 2 (count 2) | ..."
#'
#' @param counts Numeric vector. The count values.
#' @param labels Character vector. The names corresponding to the counts.
#' @param text_color Character. Color for the label text.
#' @param num_color Character. Color for the count numbers inside the brackets.
#'
#' @export
print_inline_vec <- function(counts,
                             labels,
                             text_color = "grey50",
                             num_color = "grey90") {

  # 1. Safety Check
  if (length(counts) != length(labels)) {
    stop("Vectors 'counts' and 'labels' must be the same length.")
  }

  txt_style  <- cli::make_ansi_style(text_color)
  num_style  <- cli::make_ansi_style(num_color)
  pipe_style <- cli::make_ansi_style("white")

  formatted_counts <- format(counts, big.mark = ",", trim = TRUE)

  elements <- paste0(
    txt_style(labels),
    txt_style("("),
    num_style(formatted_counts),
    txt_style(")")
  )

  final_str <- paste(elements, collapse = pipe_style(" | "))

  cli::cli_text(final_str)
}



#' Print a styled alert message
#'
#' A flexible wrapper for \code{cli} alert functions that uses \code{cli}'s
#' theming system for consistent styling. Supports glue-style interpolation
#' including \code{cli} in-line text formatting for the following:
#' \code{{.val ...}}; \code{{.arg ...}}; \code{{.field ...}} and \code{{.path ...}}.
#' Moreover, each of these in-line values is rendered in a
#' distinct highlight color (defaulting to \code{"grey98"}),
#' making it easy to draw attention to specific variables or values within the message.
#' The alert type controls the leading icon, while the base text color and font face
#' can be customized for further emphasis.
#'
#'
#' @param text Character string. The message to display. Supports \code{cli}
#'   glue-style interpolation and inline markup (e.g. \code{{.val x}},
#'   \code{{.emph x}}, \code{{?singular/plural}}).
#' @param type Character string. Alert type controlling the leading icon:
#'   \describe{
#'     \item{\code{"s"}}{Success (green tick)}
#'     \item{\code{"d"}}{Danger (red cross)}
#'     \item{\code{"w"}}{Warning (orange exclamation)}
#'     \item{\code{"i"}}{Info (blue/cyan 'i')}
#'     \item{\code{"n"}}{Plain text, no icon (default)}
#'   }
#'   Any other value falls back to a generic alert (arrow).
#' @param face Character string. Font face applied to the base text:
#'   \describe{
#'     \item{\code{"n"}}{Normal / plain (default)}
#'     \item{\code{"b"}}{Bold}
#'     \item{\code{"i"}}{Italic}
#'     \item{\code{"bi"}}{Bold italic}
#'   }
#' @param color Character string. Base text color applied to the message body.
#'   Accepts standard color names (e.g. \code{"red"}) or hex codes (e.g.
#'   \code{"#FF5733"}). Defaults to \code{"grey80"}.
#' @param highlight_color Character string. Color applied to \code{{.val ...}};
#' \code{{.arg ...}}; \code{{.field ...}} and \code{{.path ...}} spans.
#' Each of the aforementioned spans are rendered using their default themes:
#' only the colour is changed according to this argument which accepts the same
#'   values as \code{color} and defaults to \code{"grey98"}.
#' @param .envir Environment used for glue interpolation. Defaults to the
#'   caller's environment so that variables referenced in \code{text} are
#'   resolved from the calling scope.
#'
#' @details
#' Styling is applied via \code{\link[cli]{start_app}} themes rather than
#' pre-wrapping the text in ANSI codes. This ensures that \code{cli}'s own
#' inline markup (\code{{.val}}, \code{{.emph}}, etc.) composes correctly
#' with the user-specified base color and face.
#'
#' The \code{"n"} (plain) alert type uses \code{\link[cli]{cli_text}} under
#' the \code{body} theme class, producing output with no icon and no leading
#' whitespace — visually identical to the other alert types minus the bullet
#' character.
#'
#' @return Invisible \code{NULL}; called for its side effect of printing a
#'   formatted message to the console.
#'
#' @examples
#' \dontrun{
#'   # Plain message — no icon, default grey
#'   print_alert("Starting process...")
#'
#'   # Highlight interpolated values
#'   n <- 5
#'   print_alert("Loaded {.val {n}} sample{?s}", type = "s", color = "green")
#'
#'   # Warning with italic text and custom highlight
#'   print_alert("Column {.val {col}} has missing values",
#'               type = "w", color = "orange", face = "i",
#'               highlight_color = "yellow")
#'
#'   # Bold italic danger alert
#'   print_alert("Fatal: {.val {msg}}", type = "d", color = "red", face = "bi")
#' }
#'
#' @importFrom cli start_app stop_app cli_text cli_alert cli_alert_success
#'   cli_alert_danger cli_alert_warning cli_alert_info
#' @export
print_alert <- function(text,
                        type = "n",
                        face = "n",
                        color = "grey80",
                        highlight_color = "grey98",
                        .envir = parent.frame()) {

  base_style <- list(
    color         = color,
    `font-style`  = if (face %in% c("i", "bi")) "italic",
    `font-weight` = if (face %in% c("b", "bi")) "bold"
  )

  theme <- list(
    .alert             = base_style,
    `.alert-success`   = base_style,
    `.alert-danger`    = base_style,
    `.alert-warning`   = base_style,
    `.alert-info`      = base_style,
    body               = base_style,
    span.val           = list(color = highlight_color, fmt = function(x) x),
    span.arg           = list(color = highlight_color, fmt = function(x) x),
    span.field         = list(color = highlight_color, fmt = function(x) x),
    span.path          = list(color = highlight_color, fmt = function(x) x)
  )

  app <- cli::start_app(theme = theme, .auto_close = FALSE)
  on.exit(cli::stop_app(app), add = TRUE)

  dispatch <- switch(type,
                     "s" = cli::cli_alert_success,
                     "d" = cli::cli_alert_danger,
                     "w" = cli::cli_alert_warning,
                     "i" = cli::cli_alert_info,
                     "n" = cli::cli_text,
                     cli::cli_alert
  )

  dispatch(text, .envir = .envir)
  invisible(NULL)
}



#' Format text for the CLI progress steps
#'
#' A helper function that applies ANSI color and font styles (bold/italic) to a
#' text string. This is designed to be used inside \code{cli::cli_progress_step()}
#' to style the status message while preserving the spinner's functionality.
#'
#' @param text Character string. The message text to format.
#' @param color Character string. The text color. Accepts standard names (e.g.,
#'   "red", "blue") or hex codes. Defaults to "grey60".
#' @param face Character string. A shorthand code for the font face:
#'   \itemize{
#'     \item \code{"n"}: Normal (Default)
#'     \item \code{"b"}: Bold
#'     \item \code{"i"}: Italic
#'     \item \code{"bi"}: Bold Italic
#'   }
#'
#' @return A character string containing the ANSI-styled text, ready to be
#'   passed to CLI printing functions.
#'
#' @importFrom cli make_ansi_style style_bold style_italic
#'
#' @examples
#' \dontrun{
#'   # Standard usage inside a progress step
#'   cli::cli_progress_step(step_fmt("Loading data...", color = "blue"))
#'
#'   # Bold and Red for critical steps
#'   cli::cli_progress_step(step_fmt("Critical Check", color = "red", face = "b"))
#'
#'   # Italic Grey (Subtle)
#'   cli::cli_progress_step(step_fmt("Calculating metrics", face = "i"))
#' }
#' @export
step_fmt <- function(text, color = "grey60", face = "n") {

  txt_style <- cli::make_ansi_style(color)
  styled_msg <- txt_style(text)

  switch(face,
         "b"  = cli::style_bold(styled_msg),
         "i"  = cli::style_italic(styled_msg),
         "bi" = cli::style_bold(cli::style_italic(styled_msg)),
         "n"  = styled_msg,
         styled_msg # Fallback
  )
}



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
    res <- suppressWarnings(
      suppressMessages(expr)
    )
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
