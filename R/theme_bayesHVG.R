#' A \code{ggplot2} theme for \code{bayesHVG}.
#'
#' @name theme_bayesHVG
#' @author Jack R. Leary
#' @importFrom ggplot2 theme_classic theme element_rect element_line element_blank
#' @description A publication-ready theme for creating gene expression scatterplots, embedding plots, etc.
#' @param base.size The base font size. Defaults to 12.
#' @param base.lwd The base linewidth. Defaults to 0.75.
#' @param base.family The font family to be used throughout. Defaults to "sans".
#' @param umap (Optional) If set to TRUE, removes axis text and ticks for a cleaner look. Defaults to FALSE.
#' @return A \code{ggplot2} theme.
#' @export

theme_bayesHVG <- function(base.size = 12,
                           base.lwd = 0.75,
                           base.family = "sans",
                           umap = FALSE) {
  theme_bayesHVG <- ggplot2::theme_classic(base_size = base.size,
                                           base_family = base.family,
                                           base_line_size = base.lwd,
                                           base_rect_size = base.lwd) +
                    ggplot2::theme(strip.clip = "off",
                                   strip.background = ggplot2::element_rect(linewidth = base.lwd),
                                   axis.line = ggplot2::element_line(lineend = "square"))
  if (umap) {
    theme_bayesHVG <- theme_bayesHVG +
                      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                     axis.text = ggplot2::element_blank())
  }
  return(theme_bayesHVG)
}
