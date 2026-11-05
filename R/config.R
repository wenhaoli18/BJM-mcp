# Model configuration. Edit this file to point the server at your own data / fits.
# Default: the pbc3 demo model from BJM-attribution (vignette setup).

data(pbc3, package = "BJM")

fits_path <- Sys.getenv("BJM_FITS", file.path(here, "..", "models", "fits_pbc3.rds"))
if (!file.exists(fits_path))
  stop("Fitted model not found at ", fits_path, ". Run: Rscript R/fit_pbc3.R")

CFG <- list(
  name = "pbc3",
  description = paste("Primary biliary cirrhosis (Mayo, D-penicillamine trial).",
                      "BJM with 4 biomarkers; survival sub-model Surv(years, status3) ~ age + sex."),
  data = pbc3,
  fits = readRDS(fits_path),
  id_var = "id", time_var = "year",
  event_time_var = "years", event_status_var = "status",
  baseline_vars = c("age", "sex", "drug"),
  outcome = "death or transplant (composite endpoint, status3)",
  time_unit = "years",
  markers = c("serBilir", "albumin", "prothrombin", "platelets"),
  marker_labels = c(serBilir = "serum bilirubin (log)", albumin = "albumin",
                    prothrombin = "prothrombin time (transformed)", platelets = "platelets (log)"),
  marker_scale_note = "Values are on the pbc3 dataset's transformed scales (e.g. serBilir and platelets are log-transformed).",
  surv_vars = list("Tyears1", "Tyears2", "Tyears3", "Tyears4"),
  surv_funs = list(function(x) abs(x - 1), function(x) abs(x - 3),
                   function(x) abs(x - 5), function(x) abs(x - 7)),
  default_horizon = 2,
  bandcount1 = 40, bandcount2 = 80,
  caveats = c("Research demo: no model selection or calibration check was done.",
              "Not validated for clinical decision-making.")
)
