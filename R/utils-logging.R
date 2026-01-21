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

  cli::cli_rule(left = left_label, right = right_label)

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
#' A flexible wrapper for \code{cli} alert functions. It allows you to dispatch
#' different alert types (success, danger, warning, info) using shorthand codes,
#' while simultaneously controlling text color and font face (bold/italic).
#'
#' @param text Character string. The message to display.
#' @param type Character string. The type of alert icon to display. Options:
#'   \itemize{
#'     \item \code{"s"}: Success (Green Tick)
#'     \item \code{"d"}: Danger (Red Cross)
#'     \item \code{"w"}: Warning (Orange Exclamation)
#'     \item \code{"i"}: Info (Blue/Cyan 'i') - \emph{Default}
#'     \item Any other value defaults to a generic alert (Arrow).
#'   }
#' @param face Character string. The font face styling code. Options:
#'   \itemize{
#'     \item \code{"n"}: Normal (Plain) - \emph{Default}
#'     \item \code{"b"}: Bold
#'     \item \code{"i"}: Italic
#'     \item \code{"bi"}: Bold Italic
#'  }
#' @param color Character string. The color of the text. Supports standard names
#'   (e.g., "red") or hex codes. Defaults to \code{"grey80"}.
#'
#'
#' @return No return value; prints a formatted message to the console.
#' @importFrom cli make_ansi_style cli_alert_success cli_alert_danger cli_alert_warning cli_alert_info cli_alert style_bold style_italic
#'
#' @examples
#' \dontrun{
#'   # Standard Info (Grey)
#'   print_alert("Starting process...")
#'
#'   # Success with Bold Green Text
#'   print_alert("Compilation Complete", type = "s", color = "green", face = "b")
#'
#'   # Warning with Italic Orange Text
#'   print_alert("Disk space low", type = "w", color = "orange", face = "i")
#'
#'   # Critical Error (Bold Italic Red)
#'   print_alert("Fatal Exception", type = "d", color = "red", face = "bi")
#' }
#' @export
print_alert <- function(text,
                        type = "i",
                        face = "n",
                        color = "grey80"
                        ) {

  txt_style <- cli::make_ansi_style(color)

  styled_msg <- txt_style(text)

  out_message <- switch(face,
                        "b"  = cli::style_bold(styled_msg),
                        "i"  = cli::style_italic(styled_msg),
                        "bi" = cli::style_bold(cli::style_italic(styled_msg)),
                        "n"  = styled_msg,
                        styled_msg
                        )

  switch(type,
         "s" = cli::cli_alert_success(out_message),
         "d" = cli::cli_alert_danger(out_message),
         "w" = cli::cli_alert_warning(out_message),
         "i" = cli::cli_alert_info(out_message),
         cli::cli_alert(out_message)
         )
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
