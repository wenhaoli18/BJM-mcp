# Figures returned by the MCP tools:
#   plot_evidence()  - method B evidence curves (explain_risk_change_path)
#   plot_waterfall() - method A Shapley waterfall (explain_risk_change)

suppressMessages(library(ggplot2))

CJK_FAMILIES <- c("PingFang SC", "Hiragino Sans GB", "Microsoft YaHei", "SimHei",
                  "Noto Sans CJK SC", "Noto Sans SC", "Source Han Sans SC", "Source Han Sans CN",
                  "WenQuanYi Zen Hei", "WenQuanYi Micro Hei", "Heiti SC", "STHeiti")

# First installed CJK font family, or NA. Needs a UTF-8 session to draw it.
cjk_family <- function() {
  if (!isTRUE(l10n_info()[["UTF-8"]])) return(NA_character_)
  fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e)
    sub(",.*$", "", tryCatch(suppressWarnings(system("fc-list :lang=zh family", intern = TRUE,
                                                     ignore.stderr = TRUE)), error = function(e) character(0))))
  hit <- CJK_FAMILIES[CJK_FAMILIES %in% fams]
  if (length(hit)) hit[1] else NA_character_
}

# Signed label with one decimal; avoids "-0.0".
sgn1 <- function(x) { x <- round(x, 1); x[x == 0] <- 0; sprintf("%+.1f", x) }

# Push label positions apart so neighbours are at least `gap` apart, keeping order.
spread_labels <- function(y, gap) {
  o <- order(y); z <- y[o]
  for (i in seq_along(z)[-1]) z[i] <- max(z[i], z[i - 1] + gap)
  z <- z - (mean(z) - mean(y[o]))          # recentre around the original positions
  out <- y; out[o] <- z; out
}

TXT <- list(
  en = list(
    title = "Patient %s: what the new measurements say about event timing",
    subtitle = paste0("New measurements between %s and %s %s. Above 0: makes an event at that time more likely; ",
                      "below 0: less likely.\nRisk within %s %s: %.1f%% → %.1f%% (%+.1f pp). ",
                      "Labels: each marker's contribution to that change."),
    window = "%s-%s risk window", x = "Time of event (%s)", y = "Evidence (log-likelihood ratio, centred)",
    time = "time elapsed", caption = "BJM, likelihood-ratio path decomposition (method B). Research use only."),
  zh = list(
    title = "病人 %s：新测量值对事件发生时间的证据",
    subtitle = paste0("新测量：%s 到 %s %s。曲线高于 0 表示新数据让事件更可能发生在该时间，低于 0 表示更不可能。\n",
                      "%s %s内风险：%.1f%% → %.1f%%（%+.1f 个百分点）。标签为各标志物对该变化的贡献。"),
    window = "%s %s 风险窗", x = "事件发生时间（%s）", y = "证据（对数似然比，已中心化）",
    time = "时间推进", caption = "BJM 似然比路径分解（方法 B）。仅供研究，不作为临床决策依据。")
)

plot_evidence <- function(res, file, lang = "zh", tail_len = 2) {
  fam <- if (lang == "zh") cjk_family() else NA_character_
  if (lang == "zh" && is.na(fam)) lang <- "en"
  tx <- TXT[[lang]]
  if (!isTRUE(l10n_info()[["UTF-8"]])) tx$subtitle <- gsub("\u2192", "->", tx$subtitle)
  M <- CFG$markers
  lab <- if (lang == "zh" && !is.null(CFG$marker_labels_zh)) CFG$marker_labels_zh[M] else CFG$marker_labels_short[M]
  cv <- res$evidence_curves
  p0 <- cv$posterior_no_new_data
  to <- cv$to; to[!is.finite(to)] <- max(cv$from) + tail_len
  pp <- setNames(res$contributions$psi_B_pp, res$contributions$component)
  unit <- CFG$time_unit_label[[lang]]

  e <- do.call(rbind, lapply(M, function(m) {
    y <- cv[[m]] - sum(p0 * cv[[m]])
    data.frame(key = m, name = sprintf("%s  %s pp", lab[[m]], sgn1(pp[[m]])),
               t0 = cv$from, t1 = to, y = y)
  }))
  ord <- M[order(-abs(pp[M]))]
  e$name <- factor(e$name, levels = unique(e$name[match(ord, e$key)]))
  cols <- setNames(CFG$marker_colours[M], unique(e$name[match(M, e$key)]))
  endlab <- e[e$t1 == max(e$t1), ]
  endlab$ylab <- spread_labels(endlab$y, 0.065 * diff(range(e$y)))
  s2 <- res$to_time; H <- res$horizon
  fmt <- function(x) formatC(x, format = "f", digits = 1)
  fmtH <- function(x) format(x, trim = TRUE)

  g <- ggplot(e) +
    annotate("rect", xmin = s2, xmax = s2 + H, ymin = -Inf, ymax = Inf, fill = "#efeeea") +
    annotate("text", x = s2 + H + .15, y = -Inf, hjust = 0, vjust = -.8, size = 3.1, colour = "#5f5e5a",
             label = paste("\u2190", sprintf(tx$window, fmtH(H), unit))) +
    geom_hline(yintercept = 0, colour = "#a9a8a1", linewidth = .35) +
    geom_step(aes(x = t0, y = y, colour = name), direction = "hv", linewidth = .4, alpha = .45) +
    geom_segment(aes(x = t0, xend = t1, y = y, yend = y, colour = name), linewidth = 1.15) +
    geom_text(data = endlab, aes(x = t1 + .25, y = ylab, label = name, colour = name),
              hjust = 0, size = 3.3, show.legend = FALSE) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(.01, .28))) +
    labs(title = sprintf(tx$title, res$patient_id),
         subtitle = sprintf(tx$subtitle, fmt(res$from_time), fmt(res$to_time), unit, fmtH(H), unit,
                            res$risk_from_percent, res$risk_to_percent, res$change_pp),
         x = sprintf(tx$x, unit), y = tx$y,
         caption = sprintf("%s  |  %s %s pp", tx$caption, tx$time, sgn1(pp[["time"]]))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "#4a4944", size = 9.2, lineheight = 1.15),
          plot.caption = element_text(colour = "#7a7973", size = 8),
          plot.margin = margin(12, 14, 8, 10))
  if (lang == "zh") g <- g + theme(text = element_text(family = fam))
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file, g, width = 8, height = 4.8, dpi = 160, bg = "white", device = dev)
  invisible(lang)
}

TXT_WF <- list(
  en = list(
    title = "Patient %s: why did the %s-%s risk change?",
    subtitle = paste0("Risk of %s within %s %s: %.1f%% at %s %s → %.1f%% at %s %s (%+.1f pp).\n",
                      "Each bar: how much that item moved the risk (pp). They add up exactly to the change."),
    from = "At %s %s", to = "At %s %s", time = "Time elapsed",
    up = "Raises risk", down = "Lowers risk", total = "Predicted risk", tm = "Time elapsed",
    y = "%s-%s risk (%%)",
    caption = "BJM, Shapley attribution over new measurements (method A). Research use only."),
  zh = list(
    title = "病人 %s：%s %s风险为什么变了？",
    subtitle = paste0("%s %s内%s的风险：第 %s %s %.1f%% → 第 %s %s %.1f%%（%+.1f 个百分点）。\n",
                      "每根柱子表示该项使风险变化了多少个百分点，各项相加正好等于总变化。"),
    from = "第 %s %s", to = "第 %s %s", time = "时间推进",
    up = "升高风险", down = "降低风险", total = "预测风险", tm = "时间推进",
    y = "%s %s风险（%%）",
    caption = "BJM Shapley 归因（方法 A）。仅供研究，不作为临床决策依据。")
)

plot_waterfall <- function(res, file, lang = "zh") {
  fam <- if (lang == "zh") cjk_family() else NA_character_
  if (lang == "zh" && is.na(fam)) lang <- "en"
  tx <- TXT_WF[[lang]]
  if (!isTRUE(l10n_info()[["UTF-8"]])) tx$subtitle <- gsub("→", "->", tx$subtitle)
  M <- CFG$markers
  lab <- if (lang == "zh" && !is.null(CFG$marker_labels_zh)) CFG$marker_labels_zh[M] else CFG$marker_labels_short[M]
  lab <- c(time = tx$time, lab)
  unit <- CFG$time_unit_label[[lang]]
  unit1 <- if (lang == "en" && !is.null(CFG$time_unit_singular)) CFG$time_unit_singular else unit
  fmt <- function(x) formatC(x, format = "f", digits = 1)
  fmtH <- function(x) format(x, trim = TRUE)
  H <- res$horizon
  outcome <- if (lang == "zh" && !is.null(CFG$outcome_zh)) CFG$outcome_zh else CFG$outcome_short

  cb <- res$contributions
  v <- setNames(100 * cb$contribution, cb$component)
  v <- c(v["time"], v[M][order(-abs(v[M]))])
  R1 <- 100 * res$risk_from; R2 <- 100 * res$risk_to
  d <- data.frame(name = unname(lab[names(v)]), val = unname(v),
                  kind = ifelse(names(v) == "time", "time", ifelse(v > 0, "up", "down")))
  d$end <- R1 + cumsum(d$val); d$start <- d$end - d$val
  d <- rbind(data.frame(name = sprintf(tx$from, fmt(res$from_time), unit), val = R1, kind = "total", end = R1, start = 0),
             d,
             data.frame(name = sprintf(tx$to, fmt(res$to_time), unit), val = R2, kind = "total", end = R2, start = 0))
  d$x <- seq_len(nrow(d))
  d$lbl <- ifelse(d$kind == "total", sprintf("%.1f%%", d$val), sgn1(d$val))
  ytop <- max(d$start, d$end)

  g <- ggplot(d) +
    geom_rect(aes(xmin = x - .38, xmax = x + .38, ymin = start, ymax = end, fill = kind)) +
    geom_segment(data = d[d$x < max(d$x), ], aes(x = x + .38, xend = x + .62, y = end, yend = end),
                 colour = "#b5b4ad", linewidth = .3) +
    geom_text(aes(x = x, y = pmax(start, end), label = lbl), vjust = -.55, size = 3.4, colour = "#2b2b29") +
    scale_x_continuous(breaks = d$x, labels = d$name, expand = expansion(add = .3)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, NA)) +
    scale_fill_manual(values = c(total = "#c3c2b7", time = "#8a8985", up = "#e34948", down = "#2a78d6"),
                      labels = c(total = tx$total, time = tx$tm, up = tx$up, down = tx$down),
                      breaks = c("total", "time", "up", "down"), name = NULL) +
    labs(title = sprintf(tx$title, res$patient_id, fmtH(H), unit1),
         subtitle = if (lang == "zh")
           sprintf(tx$subtitle, fmtH(H), unit, outcome, fmt(res$from_time), unit, R1,
                   fmt(res$to_time), unit, R2, res$change_pp)
         else sprintf(tx$subtitle, outcome, fmtH(H), unit, R1, fmt(res$from_time), unit,
                      R2, fmt(res$to_time), unit, res$change_pp),
         x = NULL, y = sprintf(tx$y, fmtH(H), unit1), caption = tx$caption) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "#ecebe6", linewidth = .3),
          axis.text.x = element_text(colour = "#3d3c38", size = 9.5),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "#4a4944", size = 9.2, lineheight = 1.15),
          plot.caption = element_text(colour = "#7a7973", size = 8),
          plot.margin = margin(12, 14, 8, 10))
  if (lang == "zh") g <- g + theme(text = element_text(family = fam))
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file, g, width = 8, height = 4.8, dpi = 160, bg = "white", device = dev)
  invisible(lang)
}
