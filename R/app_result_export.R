# Exports reuse prepared plots and fitted summaries without fitting models.
export_mediator_label <- function(x) {
  labels <- sub("^(PC[0-9]+)_R$", "\\1", as.character(x))
  labels[labels %in% c("all", "all_mediators")] <- "All mediators"
  labels
}

export_interval_type <- function(framework) {
  if (framework == "bayesian") "95% credible interval" else "95% bootstrap confidence interval"
}

export_join_ci <- function(rows, ci, keys) {
  if (is.null(ci) || !nrow(ci)) return(rows)
  cols <- grep("_(lwr|upr)$", names(ci), value = TRUE)
  dplyr::left_join(rows, ci[, unique(c(keys, cols)), drop = FALSE], by = keys)
}

export_normalize <- function(rows, framework, kind = c("effects", "paths"), sensitivity = FALSE) {
  kind <- match.arg(kind)
  value <- function(name, default = NA_real_) {
    if (name %in% names(rows)) rows[[name]] else rep(default, nrow(rows))
  }
  bayes <- framework == "bayesian"
  if (kind == "effects") {
    columns <- if (bayes) {
      list("Total indirect effect" = c("total_IIE_est", "total_IIE_lwr", "total_IIE_upr", "total_IIE_Pr_gt_0"),
           "Direct effect" = c("direct_est", "direct_lwr", "direct_upr", "direct_Pr_gt_0"),
           "Total effect" = c("TE_est", "TE_lwr", "TE_upr", "TE_Pr_gt_0"))
    } else if (sensitivity) {
      list("Total indirect effect" = c("NIE_adjusted", "adjusted_indirect_lwr", "adjusted_indirect_upr"),
           "Direct effect" = c("NDE_adjusted", "adjusted_direct_lwr", "adjusted_direct_upr"),
           "Total effect" = c("TE", "adjusted_total_lwr", "adjusted_total_upr"))
    } else {
      list("Total indirect effect" = c("NIE", "NIE_lwr", "NIE_upr"),
           "Direct effect" = c("NDE", "NDE_lwr", "NDE_upr"),
           "Total effect" = c("TE", "TE_lwr", "TE_upr"))
    }
  } else {
    columns <- if (bayes) {
      list(a = c("alpha_est", "alpha_lwr", "alpha_upr"),
           b = c("beta_est", "beta_lwr", "beta_upr"),
           "Indirect effect" = c("IIE_est", "IIE_lwr", "IIE_upr", "IIE_Pr_gt_0"))
    } else {
      list(a = c("alpha_X_to_M", "alpha_lwr", "alpha_upr"),
           b = c(if (sensitivity) "beta_adjusted" else "beta_M_to_Y", "beta_lwr", "beta_upr"),
           "Indirect effect" = c(if (sensitivity) "indirect_adjusted" else "indirect_component", "indirect_lwr", "indirect_upr"))
    }
  }
  dplyr::bind_rows(lapply(names(columns), function(quantity) {
    cols <- columns[[quantity]]
    lower <- value(cols[2])
    upper <- value(cols[3])
    data.frame(
      contrast = value("contrast", ""),
      scenario = if (sensitivity) value("scenario", "common") else "baseline",
      varied_mediator = if (sensitivity) {
        if (bayes) value("mediator_varied", "all") else value("varied_mediator", "all")
      } else "all",
      rho = if (sensitivity) value("rho") else NA_real_,
      mediator = if (kind == "paths") export_mediator_label(value("mediator", "")) else "All mediators",
      quantity = quantity, estimate = value(cols[1]), lower = lower, upper = upper,
      interval_excludes_zero = ifelse(is.finite(lower) & is.finite(upper), lower > 0 | upper < 0, NA),
      posterior_probability_gt_zero = if (length(cols) == 4) value(cols[4]) else NA_real_,
      feasible_fraction = if (sensitivity) {
        if (bayes) value("valid_draw_fraction") else as.numeric(value("valid_parameter", TRUE))
      } else 1,
      interval_type = export_interval_type(framework), stringsAsFactors = FALSE
    )
  })) |>
    dplyr::mutate(varied_mediator = export_mediator_label(varied_mediator))
}

export_zero_crossings <- function(curves) {
  primary <- curves |>
    dplyr::filter((scenario == "common" & quantity == "Total indirect effect") |
                    (scenario == "one_at_a_time" & quantity == "Indirect effect" & mediator == varied_mediator))
  primary |>
    dplyr::group_by(contrast, scenario, varied_mediator, mediator, quantity) |>
    dplyr::group_modify(function(d, key) {
      d <- d[order(d$rho), ]
      roots <- d$rho[is.finite(d$estimate) & abs(d$estimate) < 1e-10]
      if (nrow(d) > 1L) for (i in seq_len(nrow(d) - 1L)) {
        y <- d$estimate[c(i, i + 1L)]
        if (all(is.finite(y)) && prod(sign(y)) < 0) {
          roots <- c(roots, d$rho[i] - y[1] * diff(d$rho[c(i, i + 1L)]) / diff(y))
        }
      }
      roots <- sort(unique(roots))
      data.frame(rho = if (length(roots)) roots else NA_real_,
                 status = if (length(roots)) "Crossing in computed grid" else "No crossing in computed grid",
                 method = "Linear interpolation of adjacent finite grid estimates; not an exact root")
    }) |>
    dplyr::ungroup()
}

export_display_table <- function(df) {
  df <- as.data.frame(df)
  if ("scenario" %in% names(df)) {
    df$scenario <- dplyr::recode(df$scenario, common = "Common correlation", one_at_a_time = "One mediator", baseline = "Baseline")
  }
  labels <- c(contrast = "Contrast", scenario = "Scenario", varied_mediator = "Varied mediator",
              rho = "\u03c1", mediator = "Mediator", quantity = "Quantity", estimate = "Estimate",
              lower = "Lower 95%", upper = "Upper 95%", interval_excludes_zero = "Interval excludes zero",
              posterior_probability_gt_zero = "Posterior Pr(effect > 0)", feasible_fraction = "Feasible fraction",
              interval_type = "Interval type", status = "Status", method = "Method")
  matched <- intersect(names(df), names(labels))
  names(df)[match(matched, names(df))] <- unname(labels[matched])
  df
}

write_result_workbook <- function(bundle, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Excel export requires the openxlsx package.")
  tables <- list(
    "Analysis Setup" = bundle$setup,
    "Overall Effects" = bundle$effects,
    "Mediator Paths" = bundle$paths,
    "Sensitivity Settings" = bundle$settings,
    "Sensitivity Curves" = bundle$curves,
    "Selected Rho" = bundle$selected,
    "Zero Crossings" = export_zero_crossings(bundle$curves)
  )
  for (sheet in c("Overall Effects", "Mediator Paths")) {
    tables[[sheet]] <- tables[[sheet]][, setdiff(names(tables[[sheet]]), c("scenario", "varied_mediator", "rho", "feasible_fraction")), drop = FALSE]
  }
  for (sheet in names(tables)) {
    tables[[sheet]] <- tables[[sheet]][, setdiff(
      names(tables[[sheet]]),
      c("posterior_probability_gt_zero", "feasible_fraction")
    ), drop = FALSE]
  }
  wb <- openxlsx::createWorkbook()
  for (sheet in names(tables)) {
    df <- export_display_table(tables[[sheet]])
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df, keepNA = FALSE,
                        headerStyle = openxlsx::createStyle(fgFill = "#29478F", fontColour = "#FFFFFF", textDecoration = "bold"))
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(df), widths = 23)
    if (sheet %in% c("Analysis Setup", "Sensitivity Settings")) {
      openxlsx::setColWidths(wb, sheet, cols = 2, widths = 85)
      openxlsx::addStyle(wb, sheet, openxlsx::createStyle(wrapText = TRUE), rows = seq_len(nrow(df)) + 1L, cols = 2, gridExpand = TRUE)
    }
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

export_effect_label <- function(row) {
  if (!nrow(row) || !is.finite(row$estimate[1])) return("Not feasible / unavailable")
  est <- sprintf("%.4g", row$estimate[1])
  if (isTRUE(row$interval_excludes_zero[1])) est <- paste0(est, "*")
  if (all(is.finite(c(row$lower[1], row$upper[1])))) {
    paste0(est, "\n[", sprintf("%.4g", row$lower[1]), ", ", sprintf("%.4g", row$upper[1]), "]")
  } else paste0(est, "\nInterval unavailable")
}

export_dag_plot <- function(effects, paths, title, subtitle, framework, outcome, rho = NULL) {
  med <- unique(paths$mediator)
  med <- med[order(as.integer(sub("^PC", "", med)))]
  n <- length(med)
  ink <- "#243447"
  muted <- "#5f6f7f"
  blue <- "#274690"
  violet <- "#7c5cc4"
  border <- "#dce6f2"
  wrap <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")
  subtitle <- wrap(paste(subtitle, paste("Outcome:", outcome), sep = " | "), 135)
  header_height <- 86 + 18 * (length(strsplit(subtitle, "\n", fixed = TRUE)[[1]]) - 1)
  dag_height <- max(360, 280 + 40 * n)
  card_height <- max(200, 140 + 54 * ceiling(n / 2))
  card_top <- header_height + dag_height + 48
  width <- 1200
  height <- card_top + card_height + 70
  children <- list()
  add <- function(grob) children[[length(children) + 1L]] <<- grob
  # Use a fixed export canvas in UI-like pixels, independent of the browser viewport.
  box <- function(x, y, w, h, fill = "white", stroke = border, shadow = FALSE) {
    rect <- function(dx = 0, dy = 0, bg = fill, col = stroke) grid::roundrectGrob(
      x = (x + dx) / width, y = 1 - (y + dy) / height,
      width = w / width, height = h / height, just = c("left", "top"),
      r = grid::unit(0.06, "inches"), gp = grid::gpar(fill = bg, col = col, lwd = 0.7))
    if (shadow) add(rect(0, 3, "#eaf0f8", NA))
    add(rect())
  }
  text <- function(label, x, y, size = 14, colour = ink, bold = FALSE, align = "left") {
    add(grid::textGrob(label, x = x / width, y = 1 - y / height, just = c(align, "centre"),
      gp = grid::gpar(fontfamily = "sans", fontsize = size * 0.72, col = colour,
                     fontface = if (bold) "bold" else "plain", lineheight = 1.15)))
  }
  line <- function(x, y, xend, yend, colour = border, arrow = FALSE) {
    add(grid::segmentsGrob(x0 = x / width, y0 = 1 - y / height,
      x1 = xend / width, y1 = 1 - yend / height,
      arrow = if (arrow) grid::arrow(length = grid::unit(0.12, "inches"), type = "closed") else NULL,
      gp = grid::gpar(col = colour, fill = colour, lwd = if (arrow) 2.4 else 0.7, lineend = "round")))
  }
  value <- function(row, x, y, size = 20, colour = blue, align = "left") {
    label <- strsplit(export_effect_label(row), "\n", fixed = TRUE)[[1]]
    if (nchar(label[1]) > 18) size <- min(size, 12)
    text(label[1], x, y, size, colour, TRUE, align)
    if (length(label) > 1L) text(label[2], x, y + 21, 12, muted, align = align)
  }
  effect <- function(quantity) effects[effects$quantity == quantity, , drop = FALSE]
  path <- function(m, quantity) paths[paths$mediator == m & paths$quantity == quantity, , drop = FALSE]

  text(title, 28, 29, 24, bold = TRUE)
  text(subtitle, 28, 58 + (header_height - 86) / 2, 13, muted)
  box(28, header_height, 1144, dag_height, fill = "#f8fbff")
  bottom <- header_height + dag_height - 64
  # Route the paths through the open center so coefficient panels cannot obscure them.
  for (side in c("a", "b")) {
    x <- if (side == "a") c(252, 370, 490) else c(710, 830, 948)
    y <- if (side == "a") c(bottom - 24, bottom - 24, header_height + 106) else c(header_height + 106, bottom - 24, bottom - 24)
    colour <- if (side == "a") blue else violet
    add(grid::polylineGrob(x = x / width, y = 1 - y / height,
      arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed"),
      gp = grid::gpar(col = colour, fill = colour, lwd = 2.4, linejoin = "round", lineend = "round")))
  }
  line(252, bottom, 948, bottom, "#64748b", TRUE)

  node <- function(x, y, w, label, symbol, detail = NULL, accent = NULL) {
    box(x, y, w, 88, shadow = TRUE)
    if (!is.null(accent)) line(x + 9, y + 2, x + w - 9, y + 2, accent)
    text(label, x + w / 2, y + 21, 12, muted, TRUE, "centre")
    text(symbol, x + w / 2, y + 47, 23, ink, TRUE, "centre")
    if (!is.null(detail)) text(detail, x + w / 2, y + 71, 12, muted, align = "centre")
  }
  node(70, bottom - 44, 182, "EXPOSURE", "X")
  node(948, bottom - 44, 182, "OUTCOME", "Y")
  node(470, header_height + 18, 260, "MEDIATOR", "M", paste(med, collapse = ", "), blue)
  if (!is.null(rho)) {
    box(755, header_height + 44, 124, 34, "#f1edff")
    text(paste0("\u03c1 = ", format(rho)), 817, header_height + 61, 14, violet, TRUE, "centre")
  }
  for (symbol in c("a", "b")) {
    x <- if (symbol == "a") 95 else 865
    y <- header_height + 120
    colour <- if (symbol == "a") blue else violet
    box(x, y, 240, 45 + 40 * max(1, n), shadow = TRUE)
    text(symbol, x + 120, y + 22, 20, colour, TRUE, "centre")
    for (i in seq_along(med)) {
      row_y <- y + 51 + (i - 1) * 40
      text(med[i], x + 16, row_y, 13, muted, TRUE)
      value(path(med[i], symbol), x + 224, row_y, 18, colour, "right")
    }
  }
  box(480, bottom - 44, 240, 88, shadow = TRUE)
  text("c'", 600, bottom - 24, 18, muted, TRUE, "centre")
  value(effect("Direct effect"), 600, bottom + 1, 22, ink, "centre")

  text("Effect decomposition", 28, card_top - 23, 17, bold = TRUE)
  card_width <- (1144 - 24) / 3
  quantities <- c("Total indirect effect", "Direct effect", "Total effect")
  titles <- c("Indirect effect", "Direct effect", "Total effect")
  formulas <- c("a \u00d7 b", "c'", "Total indirect + c'")
  for (j in seq_len(3)) {
    x <- 28 + (j - 1) * (card_width + 12)
    colour <- c(blue, violet, ink)[j]
    box(x, card_top, card_width, card_height)
    text(titles[j], x + 16, card_top + 26, 16, bold = TRUE)
    text(formulas[j], x + card_width - 16, card_top + 26, 12, colour, TRUE, "right")
    box(x + 12, card_top + 49, card_width - 24, 75, c("#eef4ff", "#f1edff", "#f3f6fa")[j], NA)
    text(if (j == 1) "Total indirect" else titles[j], x + 24, card_top + 68, 12, muted)
    value(effect(quantities[j]), x + card_width - 24, card_top + 73, 22, colour, "right")
    if (j == 1) for (i in seq_along(med)) {
      col <- (i - 1) %% 2
      row <- (i - 1) %/% 2
      cell_x <- x + 16 + col * (card_width - 32) / 2
      cell_y <- card_top + 146 + row * 54
      text(med[i], cell_x, cell_y, 12, muted, TRUE)
      value(path(med[i], "Indirect effect"), cell_x + (card_width - 32) / 2 - 14, cell_y, 16, blue, "right")
    }
  }
  text(paste(export_interval_type(framework), "in brackets. * Interval excludes 0."), 28, height - 38, 12, muted)
  text("Overall effects include all fitted mediators.", 28, height - 19, 12, muted)
  result <- grid::gTree(children = do.call(grid::gList, children))
  attr(result, "export_height") <- height / 100
  result
}

export_curve_plot <- function(rows, title, subtitle, rho, framework) {
  ggplot2::ggplot(rows, ggplot2::aes(x = rho, y = estimate, colour = mediator, fill = mediator, group = mediator)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = rho, colour = "grey45", linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) + ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Residual correlation (\u03c1)", y = "Indirect effect",
                   colour = "Mediator", fill = "Mediator", caption = paste("Shading:", export_interval_type(framework), "where available.")) +
    ggplot2::scale_colour_manual(values = stats::setNames(rep(c("#274690", "#7c5cc4", "#2f80ed", "#be526f", "#64748b", "#3b8a9e"), length.out = length(unique(rows$mediator))), unique(rows$mediator))) +
    ggplot2::scale_fill_manual(values = stats::setNames(rep(c("#274690", "#7c5cc4", "#2f80ed", "#be526f", "#64748b", "#3b8a9e"), length.out = length(unique(rows$mediator))), unique(rows$mediator))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "#243447", size = 18),
                   plot.subtitle = ggplot2::element_text(colour = "#5f6f7f", margin = ggplot2::margin(b = 16)),
                   plot.caption = ggplot2::element_text(colour = "#5f6f7f", hjust = 0, margin = ggplot2::margin(t = 12)),
                   panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = "#e8eef5"),
                   panel.background = ggplot2::element_rect(fill = "#f8fbff", colour = "#dce6f2"),
                   legend.position = "bottom", plot.margin = ggplot2::margin(18, 22, 18, 22))
}

export_loading_plot <- function(positive, negative, pc) {
  rows <- rbind(transform(positive, direction = "Positive"),
                transform(negative, direction = "Negative"))
  rows$direction <- factor(rows$direction, levels = c("Positive", "Negative"))
  # Row keys keep labels distinct even for small feature sets with overlapping tails.
  rows$key <- factor(seq_len(nrow(rows)), levels = rev(seq_len(nrow(rows))))
  labels <- stats::setNames(rows$gene, as.character(rows$key))
  plot <- ggplot2::ggplot(rows, ggplot2::aes(x = abs(loading), y = key, fill = direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.5f", loading)), hjust = -0.1, size = 3) +
    ggplot2::facet_wrap(~direction, scales = "free", nrow = 1) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::scale_fill_manual(values = c(Positive = "#2f80ed", Negative = "#c2413a")) +
    ggplot2::labs(title = paste(pc, "feature loadings"), x = "Absolute loading", y = NULL,
                  subtitle = "Labels show signed loadings; colors match the application.") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(12, 20, 12, 12))
  attr(plot, "export_height") <- max(6.5, 1.8 + 0.24 * max(nrow(positive), nrow(negative)))
  plot
}

write_result_figures <- function(bundle, file, context_figures = list()) {
  if (!requireNamespace("zip", quietly = TRUE)) stop("Figure export requires the zip package.")
  directory <- tempfile("mediation-figures-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  s <- bundle$selection
  effects <- bundle$effects[bundle$effects$contrast == s$mediation_contrast, ]
  paths <- bundle$paths[bundle$paths$contrast == s$mediation_contrast, ]
  if (s$mediator != "all") paths <- paths[paths$mediator == s$mediator, ]
  subtitle <- paste(bundle$framework_label, s$mediation_contrast, sep = " | ")
  save <- function(name, plot) {
    height <- attr(plot, "export_height")
    if (is.null(height)) height <- 6.5
    ggplot2::ggsave(file.path(directory, name), plot, width = 12, height = height, units = "in", dpi = 300, bg = "white")
  }
  save("Mediation_Model.png", export_dag_plot(effects, paths, "Mediation Model", subtitle, bundle$framework, bundle$outcome))
  scenario_label <- if (s$scenario == "common") "Common correlation" else paste("One mediator:", s$varied_mediator)
  sub <- paste(bundle$framework_label, s$sensitivity_contrast, scenario_label, sep = " | ")
  curves <- bundle$curves[bundle$curves$contrast == s$sensitivity_contrast & bundle$curves$scenario == s$scenario, ]
  primary <- if (s$scenario == "common") {
    curves[curves$quantity == "Total indirect effect", ]
  } else curves[curves$quantity == "Indirect effect" & curves$mediator == s$varied_mediator & curves$varied_mediator == s$varied_mediator, ]
  save("Sensitivity_Curve.png", export_curve_plot(primary, if (s$scenario == "common") "Total Indirect Effect" else paste(s$varied_mediator, "Indirect Effect"), sub, s$rho, bundle$framework))
  selected <- bundle$selected[bundle$selected$contrast == s$sensitivity_contrast, ]
  effects <- selected[selected$mediator == "All mediators", ]
  paths <- selected[selected$mediator != "All mediators", ]
  if (s$scenario == "one_at_a_time") paths <- paths[paths$mediator == s$varied_mediator, ]
  save("Sensitivity_Decomposition.png", export_dag_plot(effects, paths, "Effects at Selected \u03c1",
       sub, bundle$framework, bundle$outcome, rho = s$rho))
  if (s$scenario == "one_at_a_time") {
    comparison <- curves[curves$quantity == "Indirect effect" & curves$mediator == curves$varied_mediator, ]
    comparison_subtitle <- paste(bundle$framework_label, s$sensitivity_contrast, "Each mediator varied separately", sep = " | ")
    save("Mediator_Comparison.png", export_curve_plot(comparison, "All Mediators Comparison", comparison_subtitle, s$rho, bundle$framework))
  }
  for (name in names(context_figures)) save(name, context_figures[[name]])
  zip::zipr(file, files = list.files(directory, full.names = TRUE), include_directories = FALSE)
}
