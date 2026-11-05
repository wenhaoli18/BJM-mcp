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

# Method B (BJM-attribution, psi.R / run.R): likelihood-ratio path decomposition.
#   Binned posterior of T at to_time is built from predictRisk() CDFs at a grid of
#   horizons, for each subset S of markers whose new data are added: P_S(bin).
#   Evidence curve of the new data: l_S(t) = log P_S(t) - log P_none(t);
#   per-marker evidence curve phi_k(t) = Shapley value of l_S(t).
#   psi_k = int_0^1 Cov_eps(1{T in window}, phi_k(T)) d eps, where the tempered
#   posterior is q_eps ∝ P_none * exp(eps * l_all). sum_k psi_k = data part of dR
#   exactly in the binned model. 'time' is the same as in method A.
tool_explain_path <- function(req) {
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
  if (!any(D[[CFG$time_var]] > s1)) stop("No new visit in (from_time, to_time]; nothing to attribute except time.")

  M <- CFG$markers; K <- length(M)
  hs <- sort(unique(c(CFG$path_horizons, H)))
  w  <- c(hs <= H, FALSE)                                   # bins inside the window (+ tail bin)
  dl_of <- function(keep, s) lapply(M, function(m) {
    x <- D[!(D[[CFG$time_var]] > s1 & !(m %in% keep)), ]; x[x[[CFG$time_var]] <= s, ] })
  risk_bc <- function(dl, s, h, b1, b2) {
    r <- suppressWarnings(predictRisk(
      dl, CFG$fits$l, CFG$fits$s, prediction_time = s, horizon = h,
      time_variable = CFG$time_var, survival_variable_all = CFG$surv_vars,
      survival_trans_function = CFG$surv_funs, bandcount1 = b1, bandcount2 = b2))
    unname(r$risk_prob_1[1]) }
  key  <- function(S) if (length(S)) paste(sort(S), collapse = "|") else "none"
  subs <- unlist(lapply(0:K, function(k) combn(M, k, simplify = FALSE)), recursive = FALSE)
  shap <- function(val) { out <- setNames(vector("list", K), M)
    for (p in M) { o <- setdiff(M, p); acc <- 0
      for (k in 0:(K - 1)) for (S in combn(o, k, simplify = FALSE))
        acc <- acc + factorial(k) * factorial(K - k - 1) / factorial(K) * (val[[key(c(S, p))]] - val[[key(S)]])
      out[[p]] <- acc }
    out }

  t0 <- Sys.time()
  cores <- max(1L, min(length(subs), as.integer(Sys.getenv("BJM_CORES", parallel::detectCores()))))
  run <- function(S) {
    dl <- dl_of(S, s2)
    Fm <- vapply(hs, function(h) risk_bc(dl, s2, h, CFG$path_bandcount1, CFG$path_bandcount2), numeric(1))
    list(F = Fm, RA = risk_bc(dl, s2, H, CFG$bandcount1, CFG$bandcount2)) }
  out <- if (cores > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(subs, run, mc.cores = cores) else lapply(subs, run)
  names(out) <- vapply(subs, key, "")
  if (any(vapply(out, inherits, NA, "try-error"))) stop("predictRisk failed in a worker.")

  P <- lapply(out, function(o) { m <- pmax(c(o$F[1], diff(o$F), 1 - o$F[length(o$F)]), 1e-12); m / sum(m) })
  ell <- lapply(P, function(p) log(p) - log(P$none)); L <- ell[[key(M)]]
  phi <- shap(ell)
  gx <- c(-0.9739065285,-0.8650633667,-0.6794095683,-0.4333953941,-0.1488743390,
           0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285)
  gw <- c(0.0666713443,0.1494513492,0.2190863625,0.2692667143,0.2955242247,
          0.2955242247,0.2692667143,0.2190863625,0.1494513492,0.0666713443)
  e <- (gx + 1) / 2; ew <- gw / 2
  psi <- vapply(M, function(p) sum(vapply(seq_along(e), function(i) {
    q <- P$none * exp(e[i] * L); q <- q / sum(q)
    ew[i] * (sum(q[w] * phi[[p]][w]) - sum(q[w]) * sum(q * phi[[p]])) }, numeric(1))), numeric(1))

  # Method A on the same subsets, for comparison
  VA <- lapply(out, `[[`, "RA"); shA <- unlist(shap(VA))
  R1 <- risk_bc(dl_of(character(0), s1), s1, H, CFG$bandcount1, CFG$bandcount2)
  time_eff <- VA$none - R1; R2 <- VA[[key(M)]]
  dR_data_bin <- sum(P[[key(M)]][w]) - sum(P$none[w])

  bins <- data.frame(from = s2 + c(0, hs), to = s2 + c(hs, Inf), in_window = w)
  curves <- data.frame(bins, posterior_no_new_data = P$none, posterior_all_new_data = P[[key(M)]],
                       evidence_total = L, sapply(M, function(p) phi[[p]]), check.names = FALSE)
  contrib <- data.frame(
    component = c("time", M),
    label = c("time elapsed / survived to later visit", CFG$marker_labels[M]),
    psi_B_pp = round(100 * c(time_eff, psi), 2),
    shapley_A_pp = round(100 * c(time_eff, shA[M]), 2),
    row.names = NULL)
  contrib <- contrib[order(-abs(contrib$psi_B_pp)), ]

  list(patient_id = id, horizon = H, from_time = s1, to_time = s2,
       risk_from_percent = round(100 * R1, 1), risk_to_percent = round(100 * R2, 1),
       change_pp = round(100 * (R2 - R1), 2),
       contributions = contrib,
       checks = list(
         sum_psi_minus_binned_dR_data = sum(psi) - dR_data_bin,
         binning_error_dR_data = dR_data_bin - (R2 - VA$none),
         closure_error_B = (R2 - R1) - time_eff - sum(psi),
         evidence_strength_sd = sqrt(sum(P$none * L^2) - sum(P$none * L)^2)),
       evidence_curves = curves,
       compute_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
       method = paste("Method B (likelihood-ratio path decomposition): the new data's log-likelihood-ratio",
                      "curve over event time is split into per-marker evidence curves phi_k(t) (Shapley);",
                      "psi_k integrates, over evidence strength eps in [0,1], the posterior covariance between",
                      "'event in window' and phi_k. 'time' is shared with method A. shapley_A_pp is method A",
                      "for comparison; the two agree to first order and differ at second order.",
                      "B uses a binned posterior from predictRisk() CDFs, so it carries a small binning error."),
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
      explain_risk_change_path = tool_explain_path,
      stop("Unknown tool: ", req$tool))
    list(ok = TRUE, result = fn(req))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  cat(toJSON(out, auto_unbox = TRUE, digits = 6, na = "null", null = "null", dataframe = "rows"))
}
main()
