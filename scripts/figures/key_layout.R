# Shared tree-figure key layout: the Figure 2 domain-key centring (714d8ab lineage), used by
# Figure 2, the full origin trees (Figure S3) and the accessory-domain trees (Figure S4).
# Layout only: no tips, trees, supports, colours or values are touched here.

# Width of the widest key label in panel data units. Text widths are measured in the same
# physical units as ggplot's text size, then converted with the panel's data-to-device scale
# (the saved plot is width_in inches wide).
key_label_units <- function(p, labels, legend_size, width_in) {
  grDevices::pdf(NULL)
  label_pt <- max(as.numeric(grid::convertWidth(
    grid::stringWidth(labels), "pt", valueOnly=TRUE))) *
    legend_size * ggplot2::.pt / 12
  panel_x_range <- diff(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range)
  grDevices::dev.off()
  label_pt * panel_x_range / (width_in * 72)
}

# Left x of a key whose complete visible extent (last swatch plus its widest label and the
# panel's right-side plotting clearance) is centred under [column_left, column_right].
centred_key_x <- function(column_left, column_right, n_items, legend_cols, col_spacing,
                          label_units, right_clearance) {
  key_width <- (min(legend_cols, n_items) - 1) * col_spacing +
    3.5 + label_units + right_clearance
  (column_left + column_right - key_width) / 2
}

# Taxon (clade) colour key in up to `cols` columns at the domain-key column spacing,
# centred under the same architecture column and enclosed in a thin frame.
add_framed_taxon_key <- function(p, clades, labels, colours, column_left, column_right,
                                 y_start, cols, col_spacing, row_h, legend_size, width_in,
                                 right_clearance, measure_p = p, swatch_width = 3.0,
                                 swatch_half_height = 0.5, left_x = NULL) {
  # measure_p: plot whose x range matches the saved figure (defaults to p).
  label_units <- key_label_units(measure_p, labels, legend_size, width_in)
  # Columns never closer than the widest label plus its swatch.
  text_offset <- swatch_width + 0.5
  col_spacing <- max(col_spacing, text_offset + label_units + 1.5)
  x0 <- centred_key_x(column_left, column_right, length(clades), cols, col_spacing,
                      label_units, right_clearance)
  if (!is.null(left_x)) x0 <- left_x
  n_rows <- ceiling(length(clades) / cols)
  for (idx in seq_along(clades)) {
    x <- x0 + ((idx - 1) %% cols) * col_spacing
    y <- y_start + ((idx - 1) %/% cols) * row_h
    ccol <- colours[clades[idx]]
    if (is.na(ccol)) ccol <- "grey50"
    p <- p + annotate("rect", xmin = x, xmax = x + swatch_width, ymin = y - swatch_half_height, ymax = y + swatch_half_height,
                      fill = ccol, alpha = 0.45, color = "grey40", linewidth = 0.2) +
      annotate("text", x = x + text_offset, y = y, label = labels[idx], hjust = 0, size = legend_size)
  }
  pad <- 0.6
  xmax <- x0 + (min(cols, length(clades)) - 1) * col_spacing + text_offset + label_units + pad
  p + annotate("rect", xmin = x0 - pad, xmax = xmax,
               ymin = y_start - swatch_half_height - pad, ymax = y_start + (n_rows - 1) * row_h + swatch_half_height + pad,
               fill = NA, color = "grey40", linewidth = 0.25)
}
