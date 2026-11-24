# =============================================================================
# SENSITIVITY ANALYSIS v2 — menarche-origin flexible parametric model (FPM)
# UK Biobank breast cancer PRS study (Application 362917)
# =============================================================================
# Design (final, per user decisions Q1-Q8):
#   - Time origin   : age at menarche (menarche-origin time scale)
#   - Delayed entry : age at UKB recruitment (left truncation)
#   - Data          : df_recode_clean.csv (main-analysis cleaned data; prevalent
#                     BC already excluded; status/end_age already computed)
#   - Covariates    : ALL main-analysis clinical-reproductive-familial factors
#                     PLUS lifestyle: BMI (21001-0.0), smoking (20116-0.0),
#                     alcohol (1558-0.0)
#   - baseline age  : NOT included (per Q8: not applicable on menarche scale)
#   - First birth   : TRUE time-varying covariate; follow-up split at age at
#                     first birth when the birth occurs during follow-up.
#                     NOTE (verified from df_recode_clean.csv): in this cohort
#                     NO first birth occurs after recruitment (0 women), so the
#                     time-varying specification is numerically equivalent to a
#                     fixed baseline indicator here; the specification is kept
#                     because it is the methodologically correct statement.
#   - Model         : rstpm2::stpm2, PH link, df = 5, Surv(tstart, tstop, event)
#                       tstart = recruitment age - menarche age
#                       tstop  = event/censor age - menarche age
#   - Reference profile: cohort-average rule (consistent with main analysis);
#                     categorical covariates set to their cohort MODE
#   - PRS representatives: main-analysis convention (median of 20-80 group;
#                     q80/q90/q95 boundaries for top groups)
# =============================================================================
# SETTINGS TO CONFIRM BEFORE RUNNING (see Section 0):
#   pna_as_category     : lifestyle "Prefer not to answer" handling
#   clinical_na_as_zero : how clinical NA (hysterectomy 28704, oophorectomy
#                         3430, OC 595, HRT 759) are handled — MUST mirror the
#                         main analysis
#   use_late_afb        : include time-varying late first birth (AFB>=30)
# =============================================================================

library(survival)
library(dplyr)
library(tidyr)
library(tibble)
library(rstpm2)
library(ggplot2)


# -----------------------------------------------------------------------------
#  cleaned data and map column names
# -----------------------------------------------------------------------------


names(df_raw)[names(df_raw) == "X2734.0.0"] <- "live_births"
names(df_raw)[names(df_raw) == "X2754.0.0"] <- "afb_age"
names(df_raw)[names(df_raw) == "21001-0.0"] <- "bmi"
names(df_raw)[names(df_raw) == "20116-0.0"] <- "smoking_status"
names(df_raw)[names(df_raw) == "1558-0.0"] <- "alcohol_freq"

needed <- c("eid", "baseline_age", "end_age", "status", "menarche_age_pref",
            "PRS", "FH", pc_vars, "live_births", "afb_age",
            "OC_ever", "HRT_ever", "hysterectomy_ever", "oophorectomy_bilateral",
            "benign_breast_bin", "menopause_status", "bmi", "smoking_status",
            "alcohol_freq")
missing_cols <- setdiff(needed, names(df_raw))
if (length(missing_cols)) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- df_raw %>%
  transmute(
    eid = as.character(eid),
    baseline_age = as.numeric(baseline_age),
    end_age = as.numeric(end_age),
    status = as.integer(status),
    menarche_age = as.numeric(menarche_age_pref),
    PRS = as.numeric(PRS),
    FH = as.numeric(FH),
    across(all_of(pc_vars), ~ as.numeric(.x)),
    live_births = as.numeric(live_births),
    afb_age = as.numeric(afb_age),
    OC_ever = as.numeric(OC_ever),
    HRT_ever = as.numeric(HRT_ever),
    hysterectomy_ever = as.numeric(hysterectomy_ever),
    oophorectomy_bilateral = as.numeric(oophorectomy_bilateral),
    benign_breast_bin = as.numeric(benign_breast_bin),
    menopause_status = as.character(menopause_status),
    bmi = as.numeric(bmi),
    smoking_status = as.numeric(smoking_status),
    alcohol_freq = as.numeric(alcohol_freq)
  )

stopifnot(nrow(df) == nrow(df_raw), !anyDuplicated(df$eid))

cat("=== Cohort QC ===\n")
print(tibble(N = nrow(df),
             events = sum(df$status == 1, na.rm = TRUE),
             baseline_age_min = min(df$baseline_age, na.rm = TRUE),
             baseline_age_max = max(df$baseline_age, na.rm = TRUE),
             t_entry_min = min(df$baseline_age - df$menarche_age, na.rm = TRUE),
             t_entry_max = max(df$baseline_age - df$menarche_age, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 2. Missing / "Prefer not to answer" counts for lifestyle variables
#    (informative; decision is set via pna_as_category in Section 0)
# -----------------------------------------------------------------------------
smoking_labels <- c("0"="Never", "1"="Previous", "2"="Current",
                    "-3"="Prefer not to answer")
alcohol_labels <- c("1"="Daily or almost daily",
                    "2"="Three or four times a week",
                    "3"="Once or twice a week",
                    "4"="One to three times a month",
                    "5"="Special occasions only",
                    "6"="Never",
                    "-3"="Prefer not to answer")

count_cat <- function(x, labels) {
  x <- as.numeric(x)
  tab <- table(factor(x, levels = as.numeric(names(labels))), useNA = "ifany")
  tibble(category = names(tab), n = as.integer(tab),
         pct = round(100 * as.integer(tab) / sum(tab), 2))
}
cat("=== Smoking status (20116) ===\n"); print(count_cat(df$smoking_status, smoking_labels))
cat("=== Alcohol intake frequency (1558) ===\n"); print(count_cat(df$alcohol_freq, alcohol_labels))
cat("=== BMI (21001) ===\n")
print(tibble(n_total = nrow(df),
             n_NA = sum(is.na(df$bmi)),
             n_negative = sum(!is.na(df$bmi) & df$bmi < 0),
             n_valid = sum(!is.na(df$bmi) & df$bmi > 0),
             mean_valid = round(mean(df$bmi[df$bmi > 0], na.rm = TRUE), 2),
             sd_valid = round(sd(df$bmi[df$bmi > 0], na.rm = TRUE), 2),
             n_lt15 = sum(!is.na(df$bmi) & df$bmi > 0 & df$bmi < 15),
             n_gt60 = sum(!is.na(df$bmi) & df$bmi > 60)))

# -----------------------------------------------------------------------------
# 3. Clean variables
# -----------------------------------------------------------------------------
df <- df %>%
  mutate(
    # continuous sentinel codes -> NA
    bmi = case_when(is.na(bmi) ~ NA_real_, bmi < 0 ~ NA_real_, TRUE ~ bmi),
    afb_age_clean = case_when(
      is.na(afb_age) ~ NA_real_,
      afb_age <= 0  ~ NA_real_,       # includes -1/-3/-4 special codes
      TRUE ~ afb_age
    )
  )

# --- lifestyle factors ---
if (pna_as_category) {
  df <- df %>%
    mutate(
      smoking_f = factor(smoking_status, levels = c(0, 1, 2, -3),
                         labels = c("Never", "Previous", "Current",
                                    "Prefer not to answer")),
      alcohol_f = factor(alcohol_freq, levels = c(1, 2, 3, 4, 5, 6, -3),
                         labels = c("Daily or almost daily",
                                    "Three or four times a week",
                                    "Once or twice a week",
                                    "One to three times a month",
                                    "Special occasions only",
                                    "Never", "Prefer not to answer"))
    )
} else {
  df <- df %>%
    mutate(
      smoking_status = ifelse(smoking_status == -3, NA_real_, smoking_status),
      alcohol_freq   = ifelse(alcohol_freq   == -3, NA_real_, alcohol_freq),
      smoking_f = factor(smoking_status, levels = c(0, 1, 2),
                         labels = c("Never", "Previous", "Current")),
      alcohol_f = factor(alcohol_freq, levels = c(1, 2, 3, 4, 5, 6),
                         labels = c("Daily or almost daily",
                                    "Three or four times a week",
                                    "Once or twice a week",
                                    "One to three times a month",
                                    "Special occasions only",
                                    "Never"))
    )
}

# --- clinical factors: NA handling must mirror the main analysis ---
if (clinical_na_as_zero) {
  message("clinical_na_as_zero = TRUE: recoding clinical NA to 0 (never).")
  df <- df %>%
    mutate(OC_ever = ifelse(is.na(OC_ever), 0, OC_ever),
           HRT_ever = ifelse(is.na(HRT_ever), 0, HRT_ever),
           hysterectomy_ever = ifelse(is.na(hysterectomy_ever), 0, hysterectomy_ever),
           oophorectomy_bilateral = ifelse(is.na(oophorectomy_bilateral), 0, oophorectomy_bilateral))
} else {
  message("clinical_na_as_zero = FALSE: clinical NA kept missing (complete-case).")
}

# --- menopause (3 complete levels in this CSV) ---
df <- df %>%
  mutate(menopause_f = factor(menopause_status,
                              levels = c("Premenopausal", "Postmenopausal",
                                         "Uncertain")))
if (any(is.na(df$menopause_f))) {
  stop("Unexpected menopause_status values: ",
       paste(unique(df$menopause_status[is.na(df$menopause_f)]), collapse = ", "))
}

# --- parity state and menarche-origin times ---
df <- df %>%
  mutate(
    # Mirror main-analysis Parity_null coding: only live_births == 0 is
    # nulliparous; everything else (incl. PNA -3, n=163) treated as parous.
    state_entry = case_when(
      !is.na(live_births) & live_births == 0 ~ "nulliparous",
      !is.na(live_births) ~ "parous",
      !is.na(afb_age_clean) & afb_age_clean <= baseline_age ~ "parous",
      !is.na(afb_age_clean) & afb_age_clean >  baseline_age ~ "nulliparous",
      TRUE ~ NA_character_
    ),
    late_afb = ifelse(!is.na(afb_age_clean),
                      as.integer(afb_age_clean >= 30), NA_integer_),
    t_entry = baseline_age - menarche_age,
    t_exit  = end_age      - menarche_age,
    t_birth = afb_age_clean - menarche_age,
    split = (state_entry %in% "nulliparous") &
      !is.na(t_birth) & t_birth > t_entry & t_birth < t_exit
  )

cat("=== Parity / split QC ===\n")
print(df %>% summarise(
  n_live_births_PNA_treated_parous = sum(live_births == -3, na.rm = TRUE),
  n_nulliparous     = sum(state_entry == "nulliparous", na.rm = TRUE),
  n_parous          = sum(state_entry == "parous", na.rm = TRUE),
  n_unknown_parity  = sum(is.na(state_entry)),
  n_parous_afb_unknown =
    sum(state_entry == "parous" & is.na(late_afb), na.rm = TRUE),
  n_split           = sum(split, na.rm = TRUE)
))

cat("=== Time / outcome QC ===\n")
print(df %>% summarise(
  missing_baseline_age = sum(is.na(baseline_age)),
  missing_menarche     = sum(is.na(menarche_age)),
  missing_end_age      = sum(is.na(end_age)),
  missing_status       = sum(is.na(status)),
  entry_at_or_before_menarche = sum(!is.na(t_entry) & t_entry <= 0),
  exit_not_after_entry = sum(!is.na(t_entry) & !is.na(t_exit) & t_exit <= t_entry),
  event_at_or_before_recruitment =
    sum(status == 1 & end_age <= baseline_age, na.rm = TRUE)
))

# -----------------------------------------------------------------------------
# 4. Build counting-process data with TRUE time-varying first birth
# -----------------------------------------------------------------------------
fixed_cols <- c("PRS", pc_vars, "FH", "menarche_age",
                "OC_ever", "HRT_ever", "hysterectomy_ever",
                "oophorectomy_bilateral", "benign_breast_bin",
                "menopause_f", "bmi", "smoking_f", "alcohol_f")

base_fixed <- df %>%
  mutate(id = row_number()) %>%
  filter(!is.na(t_entry), !is.na(t_exit), !is.na(status),
         !is.na(state_entry)) %>%
  filter(t_entry > 0, t_exit > t_entry, status %in% c(0L, 1L))

late_afb2 <- function(x) {
  if (late_afb_unknown_as_zero) ifelse(is.na(x), 0L, x) else x
}

rows_single <- base_fixed %>%
  filter(!(split %in% TRUE)) %>%
  transmute(
    id = id, tstart = t_entry, tstop = t_exit,
    event = as.integer(status),
    AFB_tv = case_when(state_entry %in% "parous"      ~ 1L,
                       state_entry %in% "nulliparous" ~ 0L,
                       TRUE ~ NA_integer_),
    late_first_birth_tv = case_when(
      state_entry %in% "parous"      ~ late_afb2(late_afb),
      state_entry %in% "nulliparous" ~ 0L,
      TRUE ~ NA_integer_),
    across(all_of(fixed_cols))
  )

rows_pre  <- base_fixed %>%
  filter(split %in% TRUE) %>%
  transmute(
    id = id, tstart = t_entry, tstop = t_birth, event = 0L,
    AFB_tv = 0L, late_first_birth_tv = 0L,
    across(all_of(fixed_cols))
  )

rows_post <- base_fixed %>%
  filter(split %in% TRUE) %>%
  transmute(
    id = id, tstart = t_birth, tstop = t_exit,
    event = as.integer(status),
    AFB_tv = 1L, late_first_birth_tv = late_afb2(late_afb),
    across(all_of(fixed_cols))
  )

tdata_base <- bind_rows(rows_single, rows_pre, rows_post) %>%
  arrange(id, tstart)

stopifnot(all(tdata_base$tstart > 0))
stopifnot(all(tdata_base$tstop > tdata_base$tstart))
stopifnot(all(tdata_base$event %in% c(0L, 1L)))

cat("=== Interval count per participant (expect all 1 in this cohort) ===\n")
print(table(table(tdata_base$id)))

# -----------------------------------------------------------------------------
# 5. Model covariates, complete-case, FPM fitting
# -----------------------------------------------------------------------------
model_covars <- c(
  "AFB_tv",
  if (use_late_afb) "late_first_birth_tv",
  "PRS", pc_vars, "FH", "menarche_age",
  "OC_ever", "HRT_ever", "hysterectomy_ever", "oophorectomy_bilateral",
  "benign_breast_bin", "menopause_f", "bmi", "smoking_f", "alcohol_f"
)

tdata_model <- tdata_base %>%
  filter(if_all(all_of(model_covars), ~ !is.na(.x)))

cat("=== Complete-case QC ===\n")
print(tibble(
  N_before_CC      = n_distinct(tdata_base$id),
  events_before_CC = sum(tdata_base$event == 1),
  N_after_CC       = n_distinct(tdata_model$id),
  events_after_CC  = sum(tdata_model$event == 1),
  intervals_after_CC = nrow(tdata_model)
))

fpm_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~",
  paste(model_covars, collapse = " + ")
))

fit_fpm <- function(k) {
  tryCatch(
    stpm2(formula = fpm_formula, data = tdata_model, df = k,
          timeVar = "tstop", time0Var = "tstart", link.type = "PH"),
    error = function(e) {
      message("stpm2 df=", k, " failed: ", conditionMessage(e)); NULL
    }
  )
}

# df is prespecified (df = 5 per protocol/manuscript); no AIC/BIC selection.
fpm_key <- paste0("df", df_final)
fpm <- fit_fpm(df_final)
if (is.null(fpm)) stop("stpm2 ", fpm_key, " failed to fit.")
stopifnot(all(is.finite(coef(fpm))), all(is.finite(vcov(fpm))))
cat("=== Final model:", fpm_key, "(prespecified, df =", df_final, ") ===\n")
# ---- Verification against Cox (fully supports left truncation + tvc) ----
cox_formula <- as.formula(paste(
  "Surv(tstart, tstop, event) ~",
  paste(model_covars, collapse = " + "), "+ cluster(id)"
))
cox_check <- coxph(cox_formula, data = tdata_model, x = TRUE)
cox_nocluster <- coxph(as.formula(paste(
  "Surv(tstart, tstop, event) ~", paste(model_covars, collapse = " + ")
)), data = tdata_model)
cat("=== AFB_tv HR: FPM vs Cox cross-check (should be close) ===\n")
print(tibble(
  FPM_AFB_tv_HR = exp(coef(fpm)[["AFB_tv"]]),
  Cox_AFB_tv_HR = exp(coef(cox_nocluster)[["AFB_tv"]])
))
cat("=== PH diagnostic (Cox, same covariates) ===\n")
print(4(cox_nocluster))

# ---- Optional: BMI time-varying coefficient (BMI strongly violates PH) ----
# Check that the PRS association is materially unchanged when BMI's
# time-dependent effect is accommodated. tvc = list(bmi = 1) models a
# time-varying BMI coefficient (linear in log(time) on the PH scale).
fpm_tvc_bmi <- tryCatch(
  stpm2(formula = fpm_formula, data = tdata_model, df = df_final,
        timeVar = "tstop", time0Var = "tstart", link.type = "PH",
        tvc = list(bmi = 1)),
  error = function(e) {
    message("stpm2 tvc(bmi) failed: ", conditionMessage(e)); NULL
  }
)
if (!is.null(fpm_tvc_bmi)) {
  cat("=== BMI time-varying coefficient model ===\n")
  print(summary(fpm_tvc_bmi))
  tvc_b <- coef(fpm_tvc_bmi)
  if ("PRS" %in% names(tvc_b)) {
    cat("PRS HR (BMI-tvc model):", round(exp(unname(tvc_b["PRS"])), 4), "\n")
    cat("PRS HR (main df5 model):", round(exp(unname(coef(fpm)["PRS"])), 4), "\n")
  }
} else {
  message("BMI-tvc model not fitted; PH caveat for BMI should be noted in the manuscript.")
}

# -----------------------------------------------------------------------------
# 6. PRS cutpoints and representative values (main-analysis convention)
# -----------------------------------------------------------------------------
prs_cut <- quantile(df$PRS, probs = c(0.20, 0.80, 0.90, 0.95),
                    na.rm = TRUE, names = FALSE, type = 7)
print(prs_cut)

if (prs_rep_use_boundary) {
  # Main-analysis convention: median of the 20-80 group; q80/q90/q95 for
  # the top groups.
  prs_rep <- tibble(
    PRS_group = c("Reference 20-80%", "Top 20%", "Top 10%", "Top 5%"),
    PRS_value = c(
      median(df$PRS[df$PRS >= prs_cut[1] & df$PRS < prs_cut[2]], na.rm = TRUE),
      prs_cut[2], prs_cut[3], prs_cut[4]
    )
  )
} else {
  # Alternative (v1): group means
  prs_rep <- tibble(
    PRS_group = c("Reference 20-80%", "Top 20%", "Top 10%", "Top 5%"),
    PRS_value = c(
      mean(df$PRS[df$PRS >= prs_cut[1] & df$PRS < prs_cut[2]], na.rm = TRUE),
      mean(df$PRS[df$PRS >= prs_cut[2]], na.rm = TRUE),
      mean(df$PRS[df$PRS >= prs_cut[3]], na.rm = TRUE),
      mean(df$PRS[df$PRS >= prs_cut[4]], na.rm = TRUE)
    )
  )
}
print(prs_rep)
stopifnot(all(is.finite(prs_rep$PRS_value)))

# -----------------------------------------------------------------------------
# 7. m_ref and cohort-average reference profile (consistent with main analysis)
# -----------------------------------------------------------------------------
m_ref <- mean(df$menarche_age, na.rm = TRUE)   # cohort mean ~12.98 (median 13)
cat("=== m_ref (cohort mean menarche age):", round(m_ref, 4), "===\n")

mode_fun <- function(x) {
  tb <- table(x, useNA = "no")
  names(tb)[which.max(tb)]
}

ref_profile <- list(
  FH                   = mean(df$FH, na.rm = TRUE),
  menarche_age         = m_ref,
  OC_ever              = mean(df$OC_ever, na.rm = TRUE),
  HRT_ever             = mean(df$HRT_ever, na.rm = TRUE),
  hysterectomy_ever    = mean(df$hysterectomy_ever, na.rm = TRUE),
  oophorectomy_bilateral = mean(df$oophorectomy_bilateral, na.rm = TRUE),
  benign_breast_bin    = mean(df$benign_breast_bin, na.rm = TRUE),
  menopause_f          = mode_fun(df$menopause_f),
  bmi                  = mean(df$bmi, na.rm = TRUE),
  smoking_f            = mode_fun(df$smoking_f),
  alcohol_f            = mode_fun(df$alcohol_f)
)
pc_ref <- setNames(
  vapply(pc_vars, function(v) mean(df[[v]], na.rm = TRUE), numeric(1)),
  pc_vars
)
cat("=== Reference profile (cohort-average rule) ===\n")
print(ref_profile)
print(pc_ref)

# -----------------------------------------------------------------------------
# 8. Empirical support of the menarche-origin time scale
# -----------------------------------------------------------------------------
time_step <- 0.25
horizon   <- 10
t_min <- floor(min(tdata_model$tstart, na.rm = TRUE) / time_step) * time_step
t_max <- floor(max(tdata_model$tstop,  na.rm = TRUE) / time_step) * time_step
t_grid <- seq(t_min, t_max, by = time_step)

n_at_risk <- vapply(t_grid, function(tt) {
  n_distinct(tdata_model$id[tdata_model$tstart < tt & tdata_model$tstop >= tt])
}, numeric(1))

candidate_t <- t_grid[t_grid + horizon <= max(t_grid) + 1e-8]

support_qc <- bind_rows(lapply(candidate_t, function(tt) {
  idx <- which(t_grid >= tt - 1e-8 & t_grid <= tt + horizon + 1e-8)
  idx_start <- which.min(abs(t_grid - tt))
  idx_end   <- which.min(abs(t_grid - (tt + horizon)))
  events_10y <- sum(tdata_model$event == 1 &
                      tdata_model$tstop > tt & tdata_model$tstop <= tt + horizon,
                    na.rm = TRUE)
  tibble(t_start = tt, t_end = tt + horizon,
         n_at_start = n_at_risk[idx_start], n_at_end = n_at_risk[idx_end],
         min_n_at_risk_10y = min(n_at_risk[idx]),
         events_in_10y = events_10y,
         age_start = tt + m_ref, age_end = tt + horizon + m_ref)
}))

min_at_risk <- 1000
min_events_10y <- 50

support_qc <- support_qc %>%
  mutate(any_risk_set_support = min_n_at_risk_10y > 0,
         strong_support = min_n_at_risk_10y >= min_at_risk &
           events_in_10y >= min_events_10y)

print(support_qc %>% summarise(
  n_candidate_windows    = n(),
  n_with_any_risk_set    = sum(any_risk_set_support),
  n_with_strong_support  = sum(strong_support),
  earliest_candidate_age = min(age_start),
  latest_candidate_age   = max(age_start)
))

supported <- support_qc %>% filter(strong_support) %>% arrange(t_start)
if (nrow(supported) == 0) stop("No age interval meets support criteria.")
supported <- supported %>%
  mutate(run_id = cumsum(c(TRUE, diff(t_start) > time_step * 1.5)))
support_runs <- supported %>%
  group_by(run_id) %>%
  summarise(n_points = n(), min_age = min(age_start),
            max_age = max(age_start), .groups = "drop") %>%
  arrange(desc(n_points))
print(support_runs)

selected_run <- support_runs$run_id[1]
age_grid_supported <- supported %>%
  filter(run_id == selected_run) %>% pull(age_start)
support_age_min <- min(age_grid_supported)
support_age_max <- max(age_grid_supported)
cat("Support age range (selected run):",
    support_age_min, "-", support_age_max, "\n")

# -----------------------------------------------------------------------------
# 9. Prediction functions (cohort-average profile, time-varying AFB path)
# -----------------------------------------------------------------------------
extract_estimate <- function(x) {
  if (is.data.frame(x) && "Estimate" %in% names(x)) return(as.numeric(x$Estimate[1]))
  if (is.matrix(x)) return(as.numeric(x[1, 1]))
  as.numeric(x)[1]
}

get_H <- function(time, afb, late, prs) {
  if (time == 0) return(0)
  nd <- data.frame(tstart = 0, tstop = time)
  nd[["AFB_tv"]] <- as.integer(afb)
  if (use_late_afb) nd[["late_first_birth_tv"]] <- as.integer(late)
  nd[["PRS"]] <- as.numeric(prs)
  for (v in pc_vars) nd[[v]] <- pc_ref[[v]]
  nd[["FH"]]         <- ref_profile[["FH"]]
  nd[["menarche_age"]] <- ref_profile[["menarche_age"]]
  nd[["OC_ever"]]    <- ref_profile[["OC_ever"]]
  nd[["HRT_ever"]]   <- ref_profile[["HRT_ever"]]
  nd[["hysterectomy_ever"]]    <- ref_profile[["hysterectomy_ever"]]
  nd[["oophorectomy_bilateral"]] <- ref_profile[["oophorectomy_bilateral"]]
  nd[["benign_breast_bin"]]    <- ref_profile[["benign_breast_bin"]]
  nd[["menopause_f"]] <- factor(ref_profile[["menopause_f"]],
                                levels = levels(tdata_model$menopause_f))
  nd[["bmi"]] <- ref_profile[["bmi"]]
  nd[["smoking_f"]] <- factor(ref_profile[["smoking_f"]],
                              levels = levels(tdata_model$smoking_f))
  nd[["alcohol_f"]] <- factor(ref_profile[["alcohol_f"]],
                              levels = levels(tdata_model$alcohol_f))
  pred <- predict(fpm, newdata = nd, type = "cumhaz", se.fit = FALSE)
  h <- extract_estimate(pred)
  if (!is.finite(h) || h < -1e-8) stop("Bad cumhaz prediction at time ", time)
  max(h, 0)
}

cumhaz_increment_path <- function(t0, t1, birth_t, prs) {
  # birth_t: NA = never parous; else first-birth time on menarche scale
  stopifnot(length(t0) == 1L, length(t1) == 1L, t0 >= 0, t1 > t0)
  if (is.na(birth_t)) {
    return(get_H(t1, 0L, 0L, prs) - get_H(t0, 0L, 0L, prs))
  }
  if (birth_t <= t0) {
    late <- as.integer(birth_t + m_ref >= 30)
    return(get_H(t1, 1L, late, prs) - get_H(t0, 1L, late, prs))
  }
  if (birth_t >= t1) {
    return(get_H(t1, 0L, 0L, prs) - get_H(t0, 0L, 0L, prs))
  }
  late <- as.integer(birth_t + m_ref >= 30)
  (get_H(birth_t, 0L, 0L, prs) - get_H(t0, 0L, 0L, prs)) +
    (get_H(t1, 1L, late, prs) - get_H(birth_t, 1L, late, prs))
}

risk10_age <- function(age, first_birth_age, prs) {
  if (!is.finite(age) || age < m_ref) return(NA_real_)
  t0 <- age - m_ref
  t1 <- age + 10 - m_ref
  birth_t <- if (is.na(first_birth_age)) NA_real_ else first_birth_age - m_ref
  delta_H <- cumhaz_increment_path(t0, t1, birth_t, prs)
  -expm1(-max(delta_H, 0))
}

# -----------------------------------------------------------------------------
# 10. Age grid, scenarios, and 10-year risk curves
# -----------------------------------------------------------------------------
age_min_model <- max(25, m_ref + 0.25)
age_grid <- seq(from = age_min_model, to = support_age_max, by = 0.25)
if (length(age_grid) < 2) stop("Prediction grid too short.")
cat("Prediction age grid:", min(age_grid), "-", max(age_grid),
    "| extrapolated below", support_age_min, "\n")

afb_scenarios <- tibble(
  AFB_scenario    = c("No live parity", "First birth at age 30"),
  first_birth_age = c(NA_real_, 30)
)

pred_curves <- bind_rows(lapply(seq_len(nrow(afb_scenarios)), function(i) {
  bind_rows(lapply(seq_len(nrow(prs_rep)), function(j) {
    risks <- vapply(age_grid, function(a) {
      risk10_age(a, afb_scenarios$first_birth_age[i], prs_rep$PRS_value[j])
    }, numeric(1))
    tibble(age = age_grid,
           AFB_scenario = afb_scenarios$AFB_scenario[i],
           PRS_group = prs_rep$PRS_group[j],
           PRS_value = prs_rep$PRS_value[j],
           risk10 = risks)
  }))
}))

if (any(!is.finite(pred_curves$risk10))) stop("Non-finite risk predictions.")
if (any(pred_curves$risk10 < -1e-8 | pred_curves$risk10 > 1 + 1e-8)) {
  stop("Risk predictions outside [0,1].")
}

p_sensitivity <- ggplot(pred_curves,
                        aes(x = age, y = risk10, colour = PRS_group)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = risk_thresholds[1], linetype = "dashed",
             colour = "grey40") +
  geom_hline(yintercept = risk_thresholds[2], linetype = "dotted",
             colour = "grey60") +
  facet_wrap(~ AFB_scenario) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Sensitivity analysis: menarche-origin FPM",
       subtitle = "Specific 10-year breast cancer risk by PRS group (dashed: 2%, dotted: 3%)",
       x = "Chronological age (years)",
       y = "Predicted 10-year absolute risk",
       colour = "PRS group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_sensitivity)
ggsave("sensitivity_v2_curves.png", p_sensitivity,
       width = 12, height = 7.5, dpi = 300)

# -----------------------------------------------------------------------------
# 11. Threshold-attainment ages (2% and 3%)
# -----------------------------------------------------------------------------
find_crossing <- function(d, threshold, expected_step = 0.25) {
  d <- d %>% arrange(age)
  if (nrow(d) < 2 || any(!is.finite(d$age)) || any(!is.finite(d$risk10))) {
    return(tibble(crossing_age = NA_real_,
                  crossing_status = "Insufficient valid predictions"))
  }
  if (any(diff(d$age) > expected_step * 1.5)) {
    return(tibble(crossing_age = NA_real_,
                  crossing_status = "Discontinuous support range"))
  }
  above <- which(d$risk10 >= threshold)
  if (!length(above)) {
    return(tibble(crossing_age = NA_real_,
                  crossing_status = "Not reached within supported age range"))
  }
  j <- above[1]
  if (j == 1L) {
    return(tibble(crossing_age = NA_real_,
                  crossing_status =
                    "Already above threshold at lower support boundary"))
  }
  age0 <- d$age[j - 1L]; age1 <- d$age[j]
  risk0 <- d$risk10[j - 1L]; risk1 <- d$risk10[j]
  age_cross <- if (risk1 == risk0) age1 else
    age0 + (threshold - risk0) * (age1 - age0) / (risk1 - risk0)
  tibble(crossing_age = age_cross,
         crossing_status = "Crossed within supported age range")
}

crossing_results <- bind_rows(lapply(risk_thresholds, function(thr) {
  pred_curves %>%
    group_by(AFB_scenario, PRS_group, PRS_value) %>%
    group_modify(~ find_crossing(.x, threshold = thr,
                                 expected_step = 0.25)) %>%
    ungroup() %>%
    mutate(threshold = thr,
           crossing_status = case_when(
             is.na(crossing_age) ~ crossing_status,
             crossing_age < support_age_min ~
               "Model-extrapolated: below empirical support",
             crossing_age <= support_age_max ~
               "Within selected empirical-support range",
             TRUE ~ "Beyond selected empirical-support range"))
}))

print(crossing_results, n = Inf)

# -----------------------------------------------------------------------------
# 12. Extract HRs for ALL model covariates (incl. factor levels)
# -----------------------------------------------------------------------------
b_all <- coef(fpm)
cat("=== All coefficient names (verify baseline terms filtered below) ===\n")
print(names(b_all))

extract_fpm_hr <- function(fit, model_name) {
  b_all <- coef(fit)
  V <- vcov(fit)
  keep <- !startsWith(names(b_all), "nsx") & names(b_all) != "(Intercept)"
  terms <- names(b_all)[keep]
  bind_rows(lapply(terms, function(term) {
    beta <- unname(b_all[term]); se <- sqrt(V[term, term])
    data.frame(Model = model_name, Covariate = term, Beta = beta, SE = se,
               HR = exp(beta),
               Lower_95_CI = exp(beta - 1.96 * se),
               Upper_95_CI = exp(beta + 1.96 * se),
               P_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE))
  }))
}

hr_all <- extract_fpm_hr(fpm, fpm_key)
print(hr_all, row.names = FALSE)

# -----------------------------------------------------------------------------
# 13. Save outputs
# -----------------------------------------------------------------------------
write.csv(hr_all, "sensitivity_v2_HRs.csv", row.names = FALSE)
write.csv(crossing_results, "sensitivity_v2_crossing_ages.csv", row.names = FALSE)
saveRDS(list(tdata_model = tdata_model, fpm = fpm,
             support_qc = support_qc, prs_rep = prs_rep,
             ref_profile = ref_profile, m_ref = m_ref),
        "sensitivity_v2_fit.rds")

cat("=== Environment ===\n")
print(R.version.string)
print(packageVersion("rstpm2"))
print(packageVersion("survival"))
------------------------------------------------------------------------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(tibble)
fit <- readRDS("sensitivity_v2_fit.rds")
fpm          <- fit$fpm
fpm_tvc_bmi  <- fit$fpm_tvc_bmi
tdata_model  <- fit$tdata_model
ref_profile  <- fit$ref_profile
pc_ref       <- fit$pc_ref
m_ref        <- fit$m_ref
hr_all       <- fit$hr_all
use_late_afb <- fit$use_late_afb
pc_vars      <- paste0("PC", 1:10)
# =============================================================================
# Figure 1: BMI time-varying log hazard ratio
# -----------------------------------------------------------------------------
# The tvc model expresses log-HR(BMI) = b_bmi + b_tvc * B(log t), where B is
# the nsx(log(tstop), df = 1) basis. We recover B(t) exactly via the model's
# own predict() (difference of linear predictors for BMI vs BMI + 1), then
# compute the pointwise standard error by the delta method using vcov().
# =============================================================================
if (!is.null(fpm_tvc_bmi)) {
  b_bmi <- unname(coef(fpm_tvc_bmi)["bmi"])
  b_tvc <- unname(coef(fpm_tvc_bmi)["as.numeric(bmi):nsx(log(tstop), df = 1)"])
  V <- vcov(fpm_tvc_bmi)
  v_b  <- V["bmi", "bmi"]
  v_bt <- V["bmi", "as.numeric(bmi):nsx(log(tstop), df = 1)"]
  v_t  <- V["as.numeric(bmi):nsx(log(tstop), df = 1)",
             "as.numeric(bmi):nsx(log(tstop), df = 1)"]
  t_hr_grid <- seq(from = max(min(tdata_model$tstop), 1),
                   to = max(tdata_model$tstop), length.out = 300)
  make_nd <- function(tt, bmi_val) {
    nd <- data.frame(tstart = 0, tstop = tt)
    nd[["AFB_tv"]] <- 0L
    if (use_late_afb) nd[["late_first_birth_tv"]] <- 0L
    nd[["PRS"]] <- 0
    for (v in pc_vars) nd[[v]] <- pc_ref[[v]]
    nd[["FH"]]          <- ref_profile[["FH"]]
    nd[["menarche_age"]] <- ref_profile[["menarche_age"]]
    nd[["OC_ever"]]     <- ref_profile[["OC_ever"]]
    nd[["HRT_ever"]]    <- ref_profile[["HRT_ever"]]
    nd[["hysterectomy_ever"]]    <- ref_profile[["hysterectomy_ever"]]
    nd[["oophorectomy_bilateral"]] <- ref_profile[["oophorectomy_bilateral"]]
    nd[["benign_breast_bin"]]    <- ref_profile[["benign_breast_bin"]]
    nd[["menopause_f"]] <- factor(ref_profile[["menopause_f"]],
                                  levels = levels(tdata_model$menopause_f))
    nd[["bmi"]] <- bmi_val
    nd[["smoking_f"]] <- factor(ref_profile[["smoking_f"]],
                                levels = levels(tdata_model$smoking_f))
    nd[["alcohol_f"]] <- factor(ref_profile[["alcohol_f"]],
                                levels = levels(tdata_model$alcohol_f))
    nd
  }
  bmi_ref <- ref_profile[["bmi"]]
  lp0 <- vapply(t_hr_grid, function(tt)
    as.numeric(predict(fpm_tvc_bmi, newdata = make_nd(tt, bmi_ref), type = "lp")),
    numeric(1))
  lp1 <- vapply(t_hr_grid, function(tt)
    as.numeric(predict(fpm_tvc_bmi, newdata = make_nd(tt, bmi_ref + 1), type = "lp")),
    numeric(1))
  log_hr <- lp1 - lp0                 # log-HR per +1 kg/m2
  B_t    <- (log_hr - b_bmi) / b_tvc  # basis value recovered exactly
  se_log_hr <- sqrt(v_b + 2 * B_t * v_bt + B_t^2 * v_t)
  cat("Range of time-varying log-HR (per +1 kg/m2):",
      round(min(log_hr), 4), "to", round(max(log_hr), 4), "\n")
  bmi_hr_df <- data.frame(
    age       = t_hr_grid + m_ref,
    log_HR    = log_hr,
    log_HR_lo = log_hr - 1.96 * se_log_hr,
    log_HR_hi = log_hr + 1.96 * se_log_hr
  )
  p_bmi <- ggplot(bmi_hr_df, aes(x = age, y = log_HR)) +
    geom_ribbon(aes(ymin = log_HR_lo, ymax = log_HR_hi),
                alpha = 0.20, fill = "steelblue") +
    geom_line(linewidth = 1.0, colour = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1)) +
    labs(
      title = "Time-varying association of BMI with breast cancer risk",
      subtitle = "Log hazard ratio per 1 kg/m2 increase in BMI (menarche-origin flexible parametric model)",
      x = "Chronological age (years)",
      y = "Log hazard ratio per 1 kg/m2"
    ) +
    theme_bw(base_size = 13)
  ggsave("supp_fig_BMI_timevarying_logHR.png", p_bmi,
         width = 8, height = 5, dpi = 300)
  print(p_bmi)
} else {
  warning("fpm_tvc_bmi not available in RDS; run the tvc block first.")
}
# =============================================================================
# Figure 2: Forest plot of covariate hazard ratios
# =============================================================================
label_map <- c(
  "AFB_tv"            = "Ever parous (time-varying)",
  "late_first_birth_tv" = "Late first birth, AFB>=30 (time-varying)",
  "PRS"               = "PRS (per SD)",
  "FH"                = "Family history of breast cancer",
  "menarche_age"      = "Age at menarche (per year)",
  "OC_ever"           = "Oral contraceptive, ever",
  "HRT_ever"          = "Hormone replacement therapy, ever",
  "hysterectomy_ever" = "Hysterectomy, ever",
  "oophorectomy_bilateral" = "Bilateral oophorectomy, ever",
  "benign_breast_bin" = "Benign breast disease, ever",
  "menopause_fPostmenopausal" = "Postmenopausal (vs premenopausal)",
  "menopause_fUncertain"      = "Menopause uncertain (vs premenopausal)",
  "bmi"               = "BMI (per 1 kg/m2)",
  "smoking_fPrevious" = "Smoking: previous (vs never)",
  "smoking_fCurrent"  = "Smoking: current (vs never)",
  "alcohol_fThree or four times a week" = "Alcohol: 3-4x/week (vs daily)",
  "alcohol_fOnce or twice a week"       = "Alcohol: 1-2x/week (vs daily)",
  "alcohol_fOne to three times a month" = "Alcohol: 1-3x/month (vs daily)",
  "alcohol_fSpecial occasions only"     = "Alcohol: special occasions (vs daily)",
  "alcohol_fNever"    = "Alcohol: never (vs daily)"
)
forest_df <- hr_all %>%
  filter(Covariate != "(Intercept)", !grepl("^PC[0-9]+$", Covariate)) %>%
  mutate(label = label_map[Covariate]) %>%
  arrange(HR) %>%
  mutate(label = factor(label, levels = unique(label)))
p_forest <- ggplot(forest_df, aes(x = HR, y = label)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = Lower_95_CI, xmax = Upper_95_CI), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.3, 1.5, 1.8),
                labels = c("0.7", "0.8", "0.9", "1", "1.1", "1.3", "1.5", "1.8")) +
  labs(
    x = "Hazard ratio (95% CI, log scale)",
    y = NULL,
    title = "Sensitivity analysis: covariate associations with incident breast cancer",
    subtitle = paste0("Flexible parametric model (menarche-origin, df = 5 ",
                      ")")
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave("supp_fig_forest_HR.png", p_forest, width = 9, height = 7.5, dpi = 300)
print(p_forest)


