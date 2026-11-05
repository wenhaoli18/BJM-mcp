# BJM MCP backend: thin wrappers over the public BJM API (predictRisk).
# Called as:  Rscript backend.R '<json request>'   -> prints one JSON object.
# Does not modify BJM; only calls exported functions.

suppressMessages({
  library(BJM)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "config.R"))   # defines CFG (data, fits, markers, transforms)

# ---------------------------------------------------------------- helpers
patient_rows <- function(id) {
  d <- CFG$data[CFG$data[[CFG$id_var]] == id, ]
  if (!nrow(d)) stop(sprintf("Patient %s not found.", id))
  d[order(d[[CFG$time_var]]), ]
}

risk_one <- function(dl, s, horizon) {
  r <- suppressWarnings(predictRisk(
    dl, CFG$fits$l, CFG$fits$s, prediction_time = s, horizon = horizon,
    time_variable = CFG$time_var, survival_variable_all = CFG$surv_vars,
    survival_trans_function = CFG$surv_funs,
    bandcount1 = CFG$bandcount1, bandcount2 = CFG$bandcount2))
  unname(r$risk_prob_1[1])
}

prep <- function(d) { d[[CFG$id_var]] <- droplevels(factor(d[[CFG$id_var]])); d }

event_warning <- function(d, s) {
  et <- d[[CFG$event_time_var]][1]; st <- d[[CFG$event_status_var]][1]
  if (!is.na(st) && st == 1 && et <= s)
    sprintf("Recorded event at %.2f, before prediction time %.2f; the prediction conditions on being event-free at %.2f, which contradicts the record.", et, s, s)
  else NULL
}

last_visit_le <- function(d, s) max(d[[CFG$time_var]][d[[CFG$time_var]] <= s])

# ---------------------------------------------------------------- tools
tool_model_info <- function(req) {
  list(dataset = CFG$name, description = CFG$description,
       markers = CFG$markers, marker_labels = CFG$marker_labels,
       time_unit = CFG$time_unit, outcome = CFG$outcome,
       n_patients = length(unique(CFG$data[[CFG$id_var]])),
       caveats = CFG$caveats)
}

tool_patient_history <- function(req) {
  d <- patient_rows(as.character(req$patient_id))
  cols <- c(CFG$time_var, CFG$markers)
  list(patient_id = as.character(req$patient_id),
       baseline = as.list(d[1, CFG$baseline_vars, drop = FALSE]),
       follow_up = list(time = d[[CFG$event_time_var]][1],
                        status = as.character(d[[CFG$event_status_var]][1])),
       marker_scale_note = CFG$marker_scale_note,
       visits = d[, cols])
}

tool_predict_risk <- function(req) {
  id <- as.character(req$patient_id); d <- patient_rows(id)
  H  <- if (is.null(req$horizon)) CFG$default_horizon else as.numeric(req$horizon)
  s  <- if (is.null(req$prediction_time)) max(d[[CFG$time_var]]) else as.numeric(req$prediction_time)
  D  <- prep(d[d[[CFG$time_var]] <= s, ])
  if (!nrow(D)) stop("No visits at or before the prediction time.")
  r  <- risk_one(rep(list(D), length(CFG$markers)), s, H)
  list(patient_id = id, prediction_time = s, horizon = H,
       risk = r, risk_percent = round(100 * r, 1),
       n_visits_used = nrow(D), last_visit_used = max(D[[CFG$time_var]]),
       interpretation = sprintf("Probability of %s within (%.2f, %.2f] %s, given event-free at %.2f and biomarker history up to %.2f.",
                                CFG$outcome, s, s + H, CFG$time_unit, s, s),
       warning = event_warning(d, s))
}

# Shapley attribution of R(s2) - R(s1):
#   time   = R(s2 | data <= s1) - R(s1 | data <= s1)
#   marker = Shapley value over subsets of markers whose (s1, s2] data are added
tool_explain_change <- function(req) {
  id <- as.character(req$patient_id); d <- patient_rows(id)
  tv <- d[[CFG$time_var]]
  H  <- if (is.null(req$horizon)) CFG$default_horizon else as.numeric(req$horizon)
  s2 <- if (is.null(req$to_time)) max(tv) else as.numeric(req$to_time)
  s1 <- if (is.null(req$from_time)) {
    prev <- tv[tv < last_visit_le(d, s2)]
    if (!length(prev)) stop("Need at least two visits to explain a change.")
    max(prev)
  } else as.numeric(req$from_time)
  if (s1 >= s2) stop("from_time must be earlier than to_time.")
  D <- prep(d[tv <= s2, ])
  if (!any(D[[CFG$time_var]] <= s1)) stop("No visit at or before from_time.")
  new_visits <- D[D[[CFG$time_var]] > s1, c(CFG$time_var, CFG$markers)]
  if (!nrow(new_visits)) stop("No new visit in (from_time, to_time]; nothing to attribute except time.")

  M <- CFG$markers; K <- length(M)
  pr <- function(keep, s) {
    dl <- lapply(M, function(m) {
      x <- D[!(D[[CFG$time_var]] > s1 & !(m %in% keep)), ]
      x[x[[CFG$time_var]] <= s, ] })
    risk_one(dl, s, H)
  }
  key <- function(S) if (length(S)) paste(sort(S), collapse = "|") else "none"
  R1 <- pr(character(0), s1)
  subs <- unlist(lapply(0:K, function(k) combn(M, k, simplify = FALSE)), recursive = FALSE)
  V <- setNames(vapply(subs, function(S) pr(S, s2), numeric(1)), vapply(subs, key, ""))
  phi <- setNames(numeric(K), M)
  for (p in M) { o <- setdiff(M, p)
    for (k in 0:(K - 1)) for (S in combn(o, k, simplify = FALSE)) {
      w <- factorial(k) * factorial(K - k - 1) / factorial(K)
      phi[p] <- phi[p] + w * (V[key(c(S, p))] - V[key(S)]) } }
  R2 <- unname(V[key(M)]); time_eff <- unname(V["none"] - R1)

  old_last <- D[D[[CFG$time_var]] <= s1, ]; old_last <- old_last[nrow(old_last), ]
  new_last <- D[nrow(D), ]
  contrib <- data.frame(
    component = c("time", M),
    label = c("time elapsed / survived to later visit", CFG$marker_labels[M]),
    contribution = c(time_eff, unname(phi)),
    contribution_pp = round(100 * c(time_eff, unname(phi)), 2),
    value_at_from = c(NA, unlist(old_last[M])),
    value_at_to = c(NA, unlist(new_last[M])),
    row.names = NULL)
  contrib <- contrib[order(-abs(contrib$contribution)), ]

  list(patient_id = id, horizon = H, from_time = s1, to_time = s2,
       risk_from = R1, risk_to = R2,
       risk_from_percent = round(100 * R1, 1), risk_to_percent = round(100 * R2, 1),
       change_pp = round(100 * (R2 - R1), 2),
       contributions = contrib,
       closure_error = R2 - R1 - time_eff - sum(phi),
       new_visits = new_visits,
       marker_scale_note = CFG$marker_scale_note,
       method = paste("Exact Shapley decomposition on the risk scale over which markers'",
                      "new measurements in (from_time, to_time] are added; 'time' is the change",
                      "from moving the prediction time forward with no new measurements.",
                      "Contributions sum exactly to the risk change."),
       warning = event_warning(d, s2))
}

# ---------------------------------------------------------------- dispatch
main <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  req <- fromJSON(if (length(a)) a[1] else paste(readLines("stdin"), collapse = ""))
  out <- tryCatch({
    fn <- switch(req$tool,
      model_info = tool_model_info, patient_history = tool_patient_history,
      predict_risk = tool_predict_risk, explain_risk_change = tool_explain_change,
      stop("Unknown tool: ", req$tool))
    list(ok = TRUE, result = fn(req))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  cat(toJSON(out, auto_unbox = TRUE, digits = 6, na = "null", null = "null", dataframe = "rows"))
}
main()
