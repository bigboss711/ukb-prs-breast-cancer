#20260917
# ---- 套件 ----
pkgs <- c("dplyr", "survival", "readr")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(dplyr); library(survival); library(readr)



# menarche valid range
min_valid_age <- 7
max_valid_age <- 30

# 分組定義：早 ≤11, 正常 12–14, 晚 ≥15
early_cut <- 11
late_cut  <- 15



# ---- 1. 檢查 eid 欄位，轉成 character ----
if (!("eid" %in% names(df_main))) stop("主資料中找不到 'eid' 欄位")
if (!("eid" %in% names(df_men)))  stop("月經資料中找不到 'eid' 欄位")
df_main$eid <- as.character(df_main$eid)
df_men$eid  <- as.character(df_men$eid)

# ---- 2. 取月經初潮優先值函數：優先 2714-0.0 -> 2714-1.0 -> 2714-2.0 -> 2714-3.0 ----
get_pref_menarche <- function(row) {
  # row: 一個 data.frame / tibble row (with columns names like '2714-0.0' ...)
  col_names <- c("2714-0.0","2714-1.0","2714-2.0","2714-3.0")
  for (cn in col_names) {
    if (!cn %in% names(row)) next
    val <- row[[cn]]
    # NA, empty string, or non-numeric: skip
    if (is.na(val)) next
    # if read as character (some UKB csv may store negatives etc as numbers or values),
    # try to coerce to numeric
    suppressWarnings(num <- as.numeric(val))
    if (is.na(num)) next
    # exclude typical UKB invalid codes: negative values (-1, -3 etc) or out of plausibility range
    if (num < min_valid_age || num > max_valid_age) next
    return(num)
  }
  return(NA_real_)
}

# ---- 3. 對月經資料逐列套用並建立 menarche_age ----
# 為效率，使用 apply over rows of df_men for selected columns
men_cols <- intersect(c("2714-0.0","2714-1.0","2714-2.0","2714-3.0"), names(df_men))
if (length(men_cols) == 0) stop("月經資料中找不到 2714-* 的欄位")

# compute menarche_age_preferred
df_men$menarche_age_pref <- apply(df_men[, men_cols, drop = FALSE], 1, function(r) {
  # r is a vector of values in the 4 columns. We mimic get_pref_menarche logic.
  for (v in r) {
    if (is.na(v) || v == "") next
    nv <- suppressWarnings(as.numeric(v))
    if (is.na(nv)) next
    if (nv < min_valid_age || nv > max_valid_age) next
    return(nv)
  }
  return(NA_real_)
})

# optional: show how many had a valid pref value
message("月經資料：有效優先 menarche_age 數量 = ", sum(!is.na(df_men$menarche_age_pref)), "/", nrow(df_men))



# ---- 6. 剔除在所有回覆都不合理者（menarche 全 NA） ----
n_before <- nrow(df_merged)
df_clean_menarche <- df_merged %>% filter(!is.na(menarche_age_pref))
n_after <- nrow(df_clean_menarche)
message("剔除 menarche 全 NA 或不合理者: 去除 ", n_before - n_after, " 筆；剩下 ", n_after, " 筆")


# ---- 7. 建立 menarche group: Early, Normal, Late ----
df_clean_menarche <- df_clean_menarche %>%
  mutate(menarche_group = case_when(
    menarche_age_pref <= early_cut ~ "Early",
    menarche_age_pref >= late_cut  ~ "Late",
    TRUE                           ~ "Normal"
  )) %>%
  mutate(menarche_group = factor(menarche_group, levels = c("Early","Normal","Late")))

----------------------------------------------------------------------------------------------------
  
   

  #Integrated master dataset
  #  ↓
  #Identify Registry + Hospital BC records
  #↓
  #Pre-baseline BC?
  #  ┌────┴────┐
  #Yes       No
  #↓         ↓
  #Exclude   Remain in cohort
  #↓
  #Post-baseline BC
  #┌────────┴────────┐
  #↓                 ↓
  #Registry C50       Hospital C50
  #40006 + 40008      41270 + 41280
  #↓                 ↓
  #Diagnosis age      Diagnosis date
  #↓
  #Convert to age
  #└────────┬────────┘
  #↓
  #earliest BC age
  #↓
  #status + end_age
  
  # ============================================================
# Breast cancer ascertainment: registry + hospital records
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# Baseline time origin
# ------------------------------------------------------------

df_clean_menarche <- df_clean_menarche %>%
  mutate(
    baseline_date = as.Date(`53-0.0`),
    baseline_age  = as.numeric(`21022-0.0`)
  )

# ============================================================
# Pre-baseline breast cancer from cancer registry
# ============================================================

registry_code_cols <- grep(
  "^40006-\\d+\\.0$",
  names(df_clean_menarche),
  value = TRUE
)

registry_age_cols <- grep(
  "^40008-\\d+\\.0$",
  names(df_clean_menarche),
  value = TRUE
)

stopifnot(
  length(registry_code_cols) == length(registry_age_cols)
)

df_clean_menarche <- df_clean_menarche %>%
  rowwise() %>%
  mutate(
    
    registry_prebaseline_bc = {
      
      codes <- c_across(
        all_of(registry_code_cols)
      )
      
      ages <- c_across(
        all_of(registry_age_cols)
      )
      
      any(
        !is.na(codes) &
          grepl("^C50", codes, ignore.case = TRUE) &
          !is.na(ages) &
          ages <= baseline_age #因为如果癌症诊断年龄与 recruitment age 相同，应当视为 baseline/prevalent disease，而不是一个 post-entry incident event。40008 本身是由 cancer diagnosis date 推导出的 cancer diagnosis age
      )
      
    }
    
  ) %>%
  ungroup()

# ============================================================
# Hospital pre-baseline BC
# ============================================================

df_clean_menarche <- df_clean_menarche %>%
  mutate(
    
    hospital_prebaseline_bc =
      !is.na(BC_date) &
      BC_date <= baseline_date
    
  )

# ============================================================
# Combined pre-baseline breast cancer#两个 source 合并判断 prevalent BC
# ============================================================

df_clean_menarche <- df_clean_menarche %>%
  mutate(
    
    prebaseline_bc =
      registry_prebaseline_bc |
      hospital_prebaseline_bc
    
  )

# ============================================================
# Exclude pre-baseline breast cancer
# ============================================================

df_exclude_prebaslineBC <- df_clean_menarche %>%
  filter(!prebaseline_bc)

# Hard check
stopifnot(
  !any(df_exclude_prebaslineBC$prebaseline_bc)
)

# ============================================================
# Participant accounting #记录到底排除了多少人
# ============================================================

n_before_prevalent_exclusion <- nrow(df_clean_menarche)

n_excluded_registry <- sum(
  df_clean_menarche$registry_prebaseline_bc,
  na.rm = TRUE
)

n_excluded_hospital <- sum(
  df_clean_menarche$hospital_prebaseline_bc,
  na.rm = TRUE
)

n_excluded_combined <- sum(
  df_clean_menarche$prebaseline_bc,
  na.rm = TRUE
)

n_after_prevalent_exclusion <- nrow(
  df_exclude_prebaslineBC
)

cat(
  "Before prevalent BC exclusion:",
  n_before_prevalent_exclusion, "\n"
)
#Before prevalent BC exclusion: 255860 
cat(
  "Registry-only/source-positive:",
  n_excluded_registry, "\n"
)
#Registry-only/source-positive: 7600 
cat(
  "Hospital-positive:",
  n_excluded_hospital, "\n"
)
#Hospital-positive: 6241
cat(
  "Unique participants excluded:",
  n_excluded_combined, "\n"
)
#Unique participants excluded: 8236 
cat(
  "Remaining participants:",
  n_after_prevalent_exclusion, "\n"
)
#Remaining participants: 247624 
stopifnot(
  n_before_prevalent_exclusion -
    n_excluded_combined ==
    n_after_prevalent_exclusion
)

# ============================================================
#  Pre-baseline BC source concordance #source-concordance table# for reference
# ============================================================

table(
  Registry = df_clean_menarche$registry_prebaseline_bc,
  Hospital = df_clean_menarche$hospital_prebaseline_bc,
  useNA = "ifany"
)

df_analysis_exclude_prebaslineBC <- df_exclude_prebaslineBC

# ============================================================
# Convert hospital diagnosis date to attained age#统一 event-time
# ============================================================
#建立 hospital diagnosis age
df_analysis_exclude_prebaslineBC <- df_analysis_exclude_prebaslineBC %>%
  mutate(
    
    hospital_bc_age =
      baseline_age +
      as.numeric(
        difftime(
          BC_date,
          baseline_date,
          units = "days"
        )
      ) / 365.25
    
  )



# ============================================================
# Registry C50: derive earliest post-baseline BC age
#    using 40006 (ICD-10) + 40008 (age at diagnosis)
# ============================================================

registry_code_cols <- grep(
  "^40006-\\d+\\.0$",
  names(df_analysis_exclude_prebaslineBC),
  value = TRUE
)

registry_age_cols <- grep(
  "^40008-\\d+\\.0$",
  names(df_analysis_exclude_prebaslineBC),
  value = TRUE
)

stopifnot(
  length(registry_code_cols) == length(registry_age_cols)
)

df_analysis_exclude_prebaslineBC <- df_analysis_exclude_prebaslineBC %>%
  rowwise() %>%
  mutate(
    registry_bc_age = {
      
      codes <- c_across(all_of(registry_code_cols))
      ages  <- c_across(all_of(registry_age_cols))
      
      idx <- which(
        !is.na(codes) &
          grepl("^C50", codes, ignore.case = TRUE) &
          !is.na(ages) &
          ages > baseline_age
      )
      
      if (length(idx) > 0) {
        min(ages[idx])
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()
# ============================================================
#  Combined incident breast cancer event age #建立 combined incident event
#mathematical :AgeBC=min(AgeBC,registry,AgeBC,hospital)
# ============================================================

df_analysis_exclude_prebaslineBC <- df_analysis_exclude_prebaslineBC %>%
  mutate(
    
    combined_bc_age = case_when(
      
      !is.na(registry_bc_age) &
        !is.na(hospital_bc_age) ~
        pmin(
          registry_bc_age,
          hospital_bc_age
        ),
      
      !is.na(registry_bc_age) ~
        registry_bc_age,
      
      !is.na(hospital_bc_age) ~
        hospital_bc_age,
      
      TRUE ~
        NA_real_
      
    )
    
  )

# ============================================================
# Recalculate follow-up and event status
# ============================================================
df_calculate_analysis_exclude_prebaslineBC <- df_analysis_exclude_prebaslineBC %>%
  mutate(
    
    # Administrative/death/loss dates:
    lost_date  = as.Date(`191-0.0`),
    
    death_date = pmin(
      as.Date(`40000-0.0`),
      as.Date(`40000-1.0`),
      na.rm = TRUE
    ),
    
    censor_date = pmin(
      lost_date,
      death_date,
      as.Date("2022-10-31"),
      na.rm = TRUE
    ),
    
    end_age =
      baseline_age +
      as.numeric(
        difftime(
          censor_date,
          baseline_date,
          units = "days"
        )
      ) / 365.25,
    
    status =
      !is.na(combined_bc_age) &
      combined_bc_age <= end_age,
    
    end_age = if_else(
      status,
      combined_bc_age,
      end_age
    ),
    
    follow_up_years =
      end_age - baseline_age
    
  )

df_calculate_analysis_exclude_prebaslineBC <- df_calculate_analysis_exclude_prebaslineBC %>%
  mutate(
    status = as.integer(status)
  )

-----------------------------------------------------------------------
  
  #summary
  df_calculate_analysis_exclude_prebaslineBC %>%
  summarise(
    N = n(),
    Registry_postbaseline =
      sum(!is.na(registry_bc_age)),
    Hospital_postbaseline =
      sum(!is.na(hospital_bc_age)),
    Combined_cases =
      sum(!is.na(combined_bc_age))
  )



-----------------------------------------------------------------------------------------------------------------
 --------------------------------------------------------------------------------------------
  
  
  df_finishedQC <- df_finishedQC %>%
  mutate(
    afb_group = case_when(
      X2754.0.0 <= 20 ~ "Early",
      X2754.0.0 >= 30 ~ "Late",
      X2754.0.0 >= 21 & X2754.0.0 <= 29 ~ "Normal",
      TRUE ~ NA_character_
    )
  )

# 設定參考組
df_finishedQC$afb_group <- factor(
  df_finishedQC$afb_group,
  levels = c("Normal", "Early", "Late")
)

table(df_finishedQC$afb_group, useNA = "ifany")

---------------------------------------------------------------------------------
  
  #找出所有違反 Cox 條件的個案
  invalid_time <- df_finishedQC %>%
  filter(
    is.na(baseline_age) |
      is.na(end_age) |
      end_age <= baseline_age
  )

nrow(invalid_time)
summary(invalid_time$end_age - invalid_time$baseline_age)

-------------------------------------------------------------------------------------
  
  # live birth 
  df_finishedQC <- df_finishedQC %>%
  mutate(
    parity = case_when(
      X2734.0.0 == 0 ~ "0",
      X2734.0.0 == 1 ~ "1",
      X2734.0.0 >= 2 ~ "2+",
      TRUE ~ NA_character_
    )
  )
# 轉成 factor 並指定參考層級
df_finishedQC$parity <- factor(df_finishedQC$parity, levels = c("0", "1", "2+"))
---------------------------------------------------------------------------------------------
  
  # ============================================================
# Family history of breast cancer
# ============================================================

fh_cols <- c(
  grep("^p20107_i[0-3]$", names(df_finishedQC), value = TRUE),
  grep("^p20110_i[0-3]$", names(df_finishedQC), value = TRUE),
  grep("^p20111_i[0-3]$", names(df_finishedQC), value = TRUE)
)

stopifnot(length(fh_cols) > 0)

fh_matrix <- sapply(
  df_finishedQC[fh_cols],
  function(x) {
    grepl(
      "Breast cancer",
      as.character(x),
      ignore.case = TRUE
    )
  }
)

df_finishedQC <- df_finishedQC %>%
  mutate(
    family_history = rowSums(
      fh_matrix,
      na.rm = TRUE
    ) > 0
  )

---------------------------------------------------------------------------------------------------
  
  df_finishedQC_ready<- df_finishedQC  

library(dplyr)
library(survival)
library(msm)        # deltamethod
library(PredictABEL) # reclassification
library(rmda)       # decision_curve, plot_decision_curve 
library(ggplot2)
library(readr)      
set.seed(12345)

# 範例 rename，請改成你資料的實際變數名稱
df_finishedQC_ready <- df_finishedQC_ready %>%
  rename(
    PRS = p26220,
    PC1 = p22009_a1,
    PC2 = p22009_a2,
    PC3 = p22009_a3,
    PC4 = p22009_a4,
    PC5 = p22009_a5,
    PC6 = p22009_a6,
    PC7 = p22009_a7,
    PC8 = p22009_a8,
    PC9 = p22009_a9,
    PC10= p22009_a10
  )

# 檢查
head(df_finishedQC_ready)
str(df_finishedQC_ready[, c("baseline_age","end_age","status","PRS","afb_group","X2734.0.0","family_history","PC1")])
-------------------------------------------------------------------------------------------      
  #1) 產生 binary exposures（兩組交互的 binary 定義）
  #註： PRS top 10%/20%/median；AFB late >=30；Parity nulliparity 0 vs ≥1）
  
  # PRS high: top 10% example
  prs_cut <- quantile(df_finishedQC_ready$PRS, probs = 0.90, na.rm = TRUE)
df_finishedQC_ready <- df_finishedQC_ready %>% mutate(
  PRS_high = ifelse(PRS >= prs_cut, 1, 0)
)

# AFB late: assume you have afb_group factor with "Early","Reference","Late"
# If you have numeric age_first_birth use: AFB_late = ifelse(age_first_birth >= 30,1,0)
if("afb_group" %in% names(df_finishedQC_ready)) {
  df_finishedQC_ready <- df_finishedQC_ready %>% mutate(AFB_late = ifelse(afb_group == "Late", 1, 0))
} else {
  # example: numeric age_first_birth
  # df_finishedQC_ready <- df_finishedQC_ready %>% mutate(AFB_late = ifelse(age_first_birth >= 30,1,0))
}

# Family history binary (assume TRUE/FALSE or 1/0)
df_finishedQC_ready <- df_finishedQC_ready %>% mutate(FH = ifelse(as.character(family_history) %in% c("TRUE","1"), 1, 0))

# Parity binary: nulliparity as exposure (1 = nulliparous)
# if X2734.0.0 values "0","1","2+" etc
if("X2734.0.0" %in% names(df_finishedQC_ready)) {
  df_finishedQC_ready <- df_finishedQC_ready %>% mutate(Parity_null = ifelse(df_finishedQC_ready$X2734.0.0 == "0", 1, 0))
} else {
  # if numeric parity:
  # df_finishedQC_ready <- df_finishedQC_ready %>% mutate(Parity_null = ifelse(parity == 0, 1, 0))
}


#check parity_null_NA
Parity_null <- df_finishedQC_ready %>%
  filter(
    is.na(Parity_null)
    
  )

nrow(Parity_null)
summary(Parity_null)

#expel
df_finishedQC_ready <- df_finishedQC_ready %>% filter(!is.na(df_finishedQC_ready$Parity_null))

df_finishedQC_ready_add <- df_finishedQC_ready

# =====================================================================
# 
# Corrected version (v2): left-join on MAIN; auto-drop index column.
# =====================================================================
library(data.table)

# ---------------------------------------------------------------------
# 1) MAIN
# ---------------------------------------------------------------------


# Drop the pandas/row-index column that fread auto-named (V1 / Unnamed: 0 / X1 ...)
drop_idx <- intersect(names(main), c("Unnamed: 0", "V1", "X1", "row.names", "index", "X"))
if (length(drop_idx) > 0) {
  main[, (drop_idx) := NULL]
  cat("Dropped index column(s):", paste(drop_idx, collapse = ", "), "\n")
} else {
  warning("No index column found. If the first column is not 'eid', stop and inspect names(main).")
}

# Safety: after dropping, eid must be the first column
stopifnot(names(main)[1] == "eid")

  #20261027
  
# ---------------------------------------------------------------------
# 0) Sentinel codes
#    -1 = Do not know;  -3 = Prefer not to answer.
#    -5 is observed in 2834/3591; meaning NOT confirmed -> treated as NA,
#    but reported separately in qc_report.
# ---------------------------------------------------------------------
sent_unk <- c(-1L, -3L, -5L)
clean_bin <- function(x) {
  x <- as.numeric(x)
  fifelse(x %in% c(0, 1), as.integer(x), NA_integer_)
}
clean_age <- function(x, lo, hi) {
  x <- as.numeric(x)
  out <- fifelse(x %in% sent_unk, NA_real_, x)
  fifelse(out >= lo & out <= hi, out, NA_real_)
}
# ---------------------------------------------------------------------
# 1) p22006  Genetic ethnic grouping (already character in this CSV)
# ---------------------------------------------------------------------
cat("\n=== 1) p22006 genetic ethnic grouping ===\n")
print(dt[, .N, by = p22006][order(-N)])
dt[, ancestry := factor(p22006)]
dt[, ancestry_known := !is.na(p22006)]
# ---------------------------------------------------------------------
# 2) 2784-0.0  Ever taken oral contraceptive pill
# ---------------------------------------------------------------------
cat("\n=== 2) 2784 ever oral contraceptive ===\n")
print(dt[, .N, by = "2784-0.0"])
dt[, OC_ever := clean_bin(colv("2784-0.0"))]
dt[, OC_ever_f := factor(OC_ever, levels = c(0L, 1L), labels = c("No", "Yes"))]
# ---------------------------------------------------------------------
# 3) X2814.0.0  Ever used HRT (instance 0)
# ---------------------------------------------------------------------
cat("\n=== 3) 2814 ever HRT ===\n")
print(dt[, .N, by = X2814.0.0])
dt[, HRT_ever := clean_bin(X2814.0.0)]
dt[, HRT_ever_f := factor(HRT_ever, levels = c(0L, 1L), labels = c("No", "Yes"))]
# ---------------------------------------------------------------------
# 4) 2724-0.0  Had menopause
#    0 = No, 1 = Yes, 2 = Not sure - had hysterectomy,
#    3 = Not sure - other reason
# ---------------------------------------------------------------------
cat("\n=== 4) 2724 had menopause ===\n")
print(dt[, .N, by = "2724-0.0"])
dt[, menopause_code := as.integer(colv("2724-0.0"))]
dt[menopause_code %in% sent_unk, menopause_code := NA_integer_]
dt[, menopause_status := factor(menopause_code,
                                levels = c(0L, 1L, 2L, 3L),
                                labels = c("Premenopausal", "Postmenopausal",
                                           "Uncertain_hysterectomy", "Uncertain_other"))]
dt[, postmenopausal := fifelse(menopause_code == 1L, 1L,
                               fifelse(menopause_code == 0L, 0L, NA_integer_))]
# ---------------------------------------------------------------------
# 5) 3581-0.0  Age at menopause (last menstrual period)
# ---------------------------------------------------------------------
cat("\n=== 5) 3581 age at menopause ===\n")
print(dt[, .N, by = "3581-0.0"][order(-N)][1:10])
dt[, age_menopause := clean_age(colv("3581-0.0"), lo = 20, hi = 80)]
dt[, age_menopause_cat := cut(age_menopause,
                              breaks = c(0, 45, 50, 55, Inf), right = FALSE,
                              labels = c("<45", "45-49", "50-54", ">=55"))]
# ---------------------------------------------------------------------
# 6) 2824-0.0  Age at hysterectomy
# ---------------------------------------------------------------------
cat("\n=== 6) 2824 age at hysterectomy ===\n")
print(dt[, .N, by = "2824-0.0"][order(-N)][1:10])
dt[, age_hysterectomy := clean_age(colv("2824-0.0"), lo = 15, hi = 80)]
# ---------------------------------------------------------------------
# 7) 2834-0.0  Bilateral oophorectomy
# ---------------------------------------------------------------------
cat("\n=== 7) 2834 bilateral oophorectomy ===\n")
print(dt[, .N, by = "2834-0.0"])
dt[, oophorectomy_bilateral := clean_bin(colv("2834-0.0"))]
dt[, oophorectomy_bilateral_f := factor(oophorectomy_bilateral,
                                        levels = c(0L, 1L), labels = c("No", "Yes"))]
# ---------------------------------------------------------------------
# 8) 3591-0.0  Ever had hysterectomy
# ---------------------------------------------------------------------
cat("\n=== 8) 3591 ever had hysterectomy ===\n")
print(dt[, .N, by = "3591-0.0"])
dt[, hysterectomy_ever := clean_bin(colv("3591-0.0"))]
dt[, hysterectomy_ever_f := factor(hysterectomy_ever,
                                   levels = c(0L, 1L), labels = c("No", "Yes"))]
# ---------------------------------------------------------------------
# 9) 3882-0.0  Age at bilateral oophorectomy
# ---------------------------------------------------------------------
cat("\n=== 9) 3882 age at bilateral oophorectomy ===\n")
print(dt[, .N, by = "3882-0.0"][order(-N)][1:10])
dt[, age_oophorectomy := clean_age(colv("3882-0.0"), lo = 15, hi = 80)]
# ---------------------------------------------------------------------
# 10) benign_breast_any_prebaseline (logical)
# ---------------------------------------------------------------------
cat("\n=== 10) benign breast disease (any pre-baseline) ===\n")
print(dt[, .N, by = benign_breast_any_prebaseline])
dt[, benign_breast_f := factor(benign_breast_any_prebaseline,
                               levels = c(FALSE, TRUE), labels = c("No", "Yes"))]
# ---------------------------------------------------------------------
# 11) Cross-checks (traceability)
# ---------------------------------------------------------------------
cat("\n=== 11) CROSS-CHECKS ===\n")
dt[, hyster_check := "OK"]
dt[hysterectomy_ever == 1L & is.na(age_hysterectomy), hyster_check := "Yes-but-no-age"]
dt[hysterectomy_ever == 0L & !is.na(age_hysterectomy), hyster_check := "No-but-has-age"]
print(dt[, .N, by = hyster_check])
dt[, ooph_check := "OK"]
dt[oophorectomy_bilateral == 1L & is.na(age_oophorectomy), ooph_check := "Yes-but-no-age"]
print(dt[, .N, by = ooph_check])
cat("2724==2 (uncertain-hysterectomy) by hysterectomy_ever:\n")
print(dt[menopause_code == 2L, .N, by = hysterectomy_ever])
dt[, benign_recompute := benign_breast_self | hospital_bdd_prebaseline]
cat("benign any vs recompute mismatches:",
    dt[, sum(benign_breast_any_prebaseline != benign_recompute, na.rm = TRUE)], "\n")
# ---------------------------------------------------------------------
# 12) Derived-variable QC
# ---------------------------------------------------------------------
cat("\n=== 12) DERIVED VARIABLES QC ===\n")
cat("PRS mean/sd/min/max:",
    dt[, mean(PRS, na.rm = TRUE)],
    dt[, sd(PRS, na.rm = TRUE)],
    dt[, min(PRS, na.rm = TRUE)],
    dt[, max(PRS, na.rm = TRUE)], "\n")
for (v in paste0("PC", 1:10)) {
  if (v %in% names(dt)) cat(v, "n_NA =", dt[, sum(is.na(get(v)))], "\n")
}
cat("family_history:\n"); print(dt[, .N, by = family_history])
cat("FH:\n");               print(dt[, .N, by = FH])
cat("parity:\n");           print(dt[, .N, by = parity])
cat("Parity_null:\n");      print(dt[, .N, by = Parity_null])
cat("AFB_late:\n");         print(dt[, .N, by = AFB_late])
cat("status:\n");           print(dt[, .N, by = status])
cat("follow_up_years range:", dt[, min(follow_up_years, na.rm = TRUE)],
    dt[, max(follow_up_years, na.rm = TRUE)], "\n")
cat("end_age < baseline_age violations:",
    dt[, sum(end_age < baseline_age, na.rm = TRUE)], "\n")
cat("prebaseline_bc:\n");   print(dt[, .N, by = prebaseline_bc])
# ---------------------------------------------------------------------
# 13) QC summary table (per raw variable)
# ---------------------------------------------------------------------
qc_summary <- function(name) {
  x <- dt[[name]]
  data.table(
    variable = name,
    n = length(x),
    n_NA = sum(is.na(x)),
    n_neg1 = sum(x == -1L, na.rm = TRUE),
    n_neg3 = sum(x == -3L, na.rm = TRUE),
    n_neg5 = sum(x == -5L, na.rm = TRUE),
    n_unique = uniqueN(x, na.rm = TRUE)
  )
}
qc_vars <- c("p22006", "2784-0.0", "X2814.0.0", "2724-0.0",
             "3581-0.0", "2824-0.0", "2834-0.0", "3591-0.0",
             "3882-0.0")
qc_report <- rbindlist(lapply(qc_vars, qc_summary), fill = TRUE)


# =====================================================================
# QC cross-checks 3591+2824 / 2834+3882  ->  PLAN A
#   Retain every "ever = Yes" woman; set unknown age to NA. No rows removed.
# Expected :
#   N = 247,618
#   hysterectomy ever=Yes = 17,794  (459 of them have unknown age)
#   oophorectomy  ever=Yes = 19,448  (498 of them have unknown age)
#   overlap = 130 | union = 827 | excluded = 0 | final N = 247,618
# =====================================================================
library(data.table)
N0 <- nrow(dt)
sent <- c(-1L, -3L, -5L)          # Do not know / Prefer not to answer / unconfirmed -5
# ---- helpers ----
to_num <- function(x) suppressWarnings(as.numeric(x))
is_yes <- function(x) { x <- to_num(x); !is.na(x) & x == 1 }
is_no  <- function(x) { x <- to_num(x); !is.na(x) & x == 0 }
valid_age <- function(x, lo, hi) {
  x <- to_num(x)
  out <- fifelse(x %in% sent, NA_real_, x)
  fifelse(out >= lo & out <= hi, out, NA_real_)
}
# =====================================================================
# Step 1: build clean EVER variables (Plan A keeps every Yes)
# =====================================================================
dt[, hysterectomy_ever := fifelse(is_yes(dt[["3591-0.0"]]), 1L,
                                  fifelse(is_no(dt[["3591-0.0"]]), 0L, NA_integer_))]
dt[, age_hysterectomy := valid_age(dt[["2824-0.0"]], lo = 15, hi = 80)]
dt[, oophorectomy_bilateral := fifelse(is_yes(dt[["2834-0.0"]]), 1L,
                                       fifelse(is_no(dt[["2834-0.0"]]), 0L, NA_integer_))]
dt[, age_oophorectomy := valid_age(dt[["3882-0.0"]], lo = 15, hi = 80)]
# =====================================================================
# Step 2: identify the cross-check flags (for reporting only; NOT exclusion)
# =====================================================================
dt[, flag_hyst_noage := (hysterectomy_ever == 1L) & is.na(age_hysterectomy)]
dt[, flag_ooph_noage := (oophorectomy_bilateral == 1L) & is.na(age_oophorectomy)]
n1 <- dt[, sum(flag_hyst_noage)]
n2 <- dt[, sum(flag_ooph_noage)]
overlap <- dt[, sum(flag_hyst_noage & flag_ooph_noage)]
union_n <- dt[, sum(flag_hyst_noage | flag_ooph_noage)]
# =====================================================================
# Step 3: REPORT (copy these numbers into the flowchart)
# =====================================================================
cat("===============================================================\n")
cat("PLAN A  (retain ever=Yes; set unknown age = NA; exclude = 0)\n")
cat("===============================================================\n")
cat("Start N                            :", N0, "\n")
cat("\n--- Cross-check 3591 + 2824 (hysterectomy) ---\n")
cat("  ever = Yes                       :", dt[, sum(hysterectomy_ever == 1L, na.rm = TRUE)], "\n")
cat("  ever = Yes, age unknown (flagged):", n1, "  -> RETAINED (age = NA)\n")
cat("\n--- Cross-check 2834 + 3882 (oophorectomy) ---\n")
cat("  ever = Yes                       :", dt[, sum(oophorectomy_bilateral == 1L, na.rm = TRUE)], "\n")
cat("  ever = Yes, age unknown (flagged):", n2, "  -> RETAINED (age = NA)\n")
cat("\n--- Overlap ---\n")
cat("  both flags true (overlap)        :", overlap, "\n")
cat("  flag1 only                       :", n1 - overlap, "\n")
cat("  flag2 only                       :", n2 - overlap, "\n")
cat("  union                            :", union_n, "\n")
cat("\n--- FINAL ---\n")
cat("  excluded                         : 0\n")
cat("  final N                          :", N0, "\n")
cat("  age_hysterectomy NA total        :", dt[, sum(is.na(age_hysterectomy))], "\n")
cat("  age_oophorectomy NA total        :", dt[, sum(is.na(age_oophorectomy))], "\n")

--------------------------------------------------------------------------------------------------------------------
  #20261029
  # =====================================================================
# Recode 7 variables per your specification (Plan: sentinels -> NA)
# Expected counts (verified from your v3 CSV):
#   2784 : No 45,833 / Yes 201,190 / NA 595
#   2814 : No 152,501 / Yes 94,358 / NA 759
#   2724 : No 59,647 / Yes 148,950 / NA 39,021
#   3581 : valid 140,155 -> early 18,315 / normal 100,557 / late 21,283; NA 107,463
#   3591 : No 201,120 / Yes 17,794 / NA 28,704
#   2834 : No 224,740 / Yes 19,448 / NA 3,430
#   benign : FALSE 239,578 / TRUE 8,040
# =====================================================================
# =====================================================================
# Recode 7 variables (v2: fixed NSE bug; uses dt[[...]] and table())
# =====================================================================
library(data.table)

N0 <- nrow(dt)

to_num <- function(x) suppressWarnings(as.numeric(x))

clean_bin <- function(x) {
  x <- to_num(x)
  fifelse(x %in% c(0, 1), as.integer(x), NA_integer_)
}
clean_age <- function(x, lo, hi) {
  x <- to_num(x)
  out <- fifelse(x %in% c(-1, -3, -5), NA_real_, x)
  fifelse(out >= lo & out <= hi, out, NA_real_)
}
report_bin <- function(name, raw, clean) {
  cat("===============================================================\n")
  cat(name, "\n")
  cat("  raw tabulation:\n")
  print(table(dt[[raw]], useNA = "ifany"))
  cat("  cleaned -> No:", sum(dt[[clean]] == 0L, na.rm = TRUE),
      "| Yes:", sum(dt[[clean]] == 1L, na.rm = TRUE),
      "| NA:", sum(is.na(dt[[clean]])), "\n")
}

# 1) 2784 Ever oral contraceptive
dt[["OC_ever"]] <- clean_bin(dt[["2784-0.0"]])
report_bin("2784 Ever oral contraceptive", "2784-0.0", "OC_ever")

# 2) 2814 Ever HRT
dt[["HRT_ever"]] <- clean_bin(dt[["X2814.0.0"]])
report_bin("2814 Ever HRT", "X2814.0.0", "HRT_ever")

# 3) 2724 Had menopause -> 0/1/NA, plus 3-level status
cat("===============================================================\n")
cat("2724 Had menopause\n")
cat("  raw tabulation:\n")
print(table(dt[["2724-0.0"]], useNA = "ifany"))
dt[["postmenopausal"]] <-
  fifelse(to_num(dt[["2724-0.0"]]) == 0, 0L,
          fifelse(to_num(dt[["2724-0.0"]]) == 1, 1L, NA_integer_))
dt[["menopause_status"]] <-
  fifelse(to_num(dt[["2724-0.0"]]) == 0, "Premenopausal",
          fifelse(to_num(dt[["2724-0.0"]]) == 1, "Postmenopausal",
                  fifelse(is.na(to_num(dt[["2724-0.0"]])), NA_character_, "Uncertain")))
cat("  postmenopausal -> No:", sum(dt[["postmenopausal"]] == 0L, na.rm = TRUE),
    "| Yes:", sum(dt[["postmenopausal"]] == 1L, na.rm = TRUE),
    "| NA:", sum(is.na(dt[["postmenopausal"]])), "\n")
print(table(dt[["menopause_status"]], useNA = "ifany"))

# 4) 3581 Age at menopause -> early/normal/late, NA retained
cat("===============================================================\n")
cat("3581 Age at menopause\n")
dt[["age_menopause"]] <- clean_age(dt[["3581-0.0"]], lo = 20, hi = 80)
dt[["menopause_age_grp"]] <- cut(
  dt[["age_menopause"]],
  breaks = c(0, 45, 55, Inf), right = FALSE,
  labels = c("Early(<45)", "Normal(45-54)", "Late(>=55)"))
print(table(dt[["menopause_age_grp"]], useNA = "ifany"))

# 5) 3591 Ever hysterectomy
dt[["hysterectomy_ever"]] <- clean_bin(dt[["3591-0.0"]])
report_bin("3591 Ever hysterectomy", "3591-0.0", "hysterectomy_ever")

# 6) 2834 Bilateral oophorectomy
dt[["oophorectomy_bilateral"]] <- clean_bin(dt[["2834-0.0"]])
report_bin("2834 Bilateral oophorectomy", "2834-0.0", "oophorectomy_bilateral")

# 7) benign_breast_any_prebaseline -> 0/1
dt[["benign_breast_bin"]] <- as.integer(dt[["benign_breast_any_prebaseline"]])
cat("===============================================================\n")
cat("benign_breast_any_prebaseline -> 0/1\n")
print(table(dt[["benign_breast_bin"]], useNA = "ifany"))



--------------------------------------------------------------------------------------------------------------------------------
  
  # =====================================================================
# 主分析 R script（final, 2026-11-17）
# 方法：絕對風險主分析 = reference profile；marginal = 補充（Supplementary Table 6）
# 所有模型使用同一組 covariates，含 hysterectomy_ever_f 三級（No/Yes/Uncertain），
# 避免 28,704 位「Not sure – had hysterectomy」被整批刪除。
# nboot 一律升到 1000。
# =====================================================================
library(survival); library(data.table); library(dplyr)
library(ggplot2); library(scales)
# 
# =====================================================================
# 1) factor / 型別準備（fcase 取代 case_when，不需 dplyr）
# =====================================================================
dt[, OC_ever_f := factor(OC_ever, levels = c(0L,1L), labels = c("No","Yes"))]
dt[, HRT_ever_f := factor(HRT_ever, levels = c(0L,1L), labels = c("No","Yes"))]
dt[, oophorectomy_f := factor(oophorectomy_bilateral, levels = c(0L,1L), labels = c("No","Yes"))]
dt[, benign_f := factor(benign_breast_bin, levels = c(0L,1L), labels = c("No","Yes"))]
dt[, menopause_status := factor(menopause_status,
                                levels = c("Premenopausal","Postmenopausal","Uncertain"))]
dt[, hysterectomy_ever_f := fcase(
  `3591-0.0` == 0, "No",
  `3591-0.0` == 1, "Yes",
  `3591-0.0` %in% c(-5,-3), "Uncertain",
  is.na(`3591-0.0`), "Uncertain",
  default = NA_character_)]
dt[, hysterectomy_ever_f := factor(hysterectomy_ever_f,
                                   levels = c("No","Yes","Uncertain"))]
dt[, PRS_high := as.integer(PRS_high)]
dt[, AFB_late  := as.integer(AFB_late)]
dt[, Parity_null := as.integer(Parity_null)]
dt[, FH := as.integer(FH)]
dt[, PRS := as.numeric(PRS)]
dt[, menarche_age_pref := as.numeric(menarche_age_pref)]
dt[, baseline_age := as.numeric(baseline_age)]
dt[, end_age := as.numeric(end_age)]
dt[, status := as.integer(status)]
dt[, (paste0("PC",1:10)) := lapply(.SD, as.numeric), .SDcols = paste0("PC",1:10)]
pc <- paste0("PC", 1:10)
# 
#=====================================================================
  # 2) complete-case（兩份：AFB 版 / Parity 版）
  # =====================================================================
cov_main <- c("baseline_age","end_age","status","PRS","AFB_late",
              "menarche_age_pref","FH",
              "OC_ever_f","HRT_ever_f","menopause_status",
              "oophorectomy_f","benign_f","hysterectomy_ever_f", pc)
analysis_afb <- na.omit(as.data.frame(dt)[, cov_main])   # parous 且 AFB 已知
cov_par <- sub("AFB_late", "Parity_null", cov_main)
analysis_par <- na.omit(as.data.frame(dt)[, cov_par])    # 全體女性
cat("AFB complete-case N =", nrow(analysis_afb), "\n")
cat("Parity complete-case N =", nrow(analysis_par), "\n")
# 
# =====================================================================
# 3) Cox：交互模型（HR 與 RERI 用）
# =====================================================================
dt_par <- dt[Parity_null == 0L & !is.na(AFB_late)]
model_prs_high_afb <- coxph(
  Surv(baseline_age, end_age, status) ~ PRS_high * AFB_late +
    menarche_age_pref + FH + OC_ever_f + HRT_ever_f +
    menopause_status + oophorectomy_f + benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = dt_par)
model_prs_afb <- coxph(
  Surv(baseline_age, end_age, status) ~ PRS * AFB_late +
    menarche_age_pref + FH + OC_ever_f + HRT_ever_f +
    menopause_status + oophorectomy_f + benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = dt_par)
model_fh_par <- coxph(
  Surv(baseline_age, end_age, status) ~ FH * Parity_null + PRS +
    menarche_age_pref + OC_ever_f + HRT_ever_f +
    menopause_status + oophorectomy_f + benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = dt)
summary(model_prs_high_afb)
summary(model_prs_afb)
summary(model_fh_par)

# 最小 vs 主模型 HR 比較（乾淨版）
extract_hr <- function(model, term = "PRS") {
  exp(cbind(HR = coef(model)[term], confint(model)[term, , drop = FALSE]))
}
model_min <- coxph(Surv(baseline_age, end_age, status) ~ PRS, data = analysis_afb)
rbind(Minimal = extract_hr(model_min), Main = extract_hr(model_prs_afb))
# VIF（car）
library(car)
vif(model_prs_afb)
vif(model_prs_high_afb)
# PH
print(cox.zph(model_prs_afb))
print(cox.zph(model_fh_par))

# =====================================================================
# 4) RERI（delta + bootstrap；定稿 nboot = 1000）
# =====================================================================
calc_RERI_from_cox <- function(cox_model, varA, varB,
                               inter_name = paste0(varA, ":", varB)) {
  if (!requireNamespace("msm", quietly = TRUE)) stop("install.packages('msm')")
  b <- coef(cox_model); V <- vcov(cox_model)
  idxA <- which(names(b) == varA); idxB <- which(names(b) == varB)
  idxI <- which(names(b) == inter_name)
  if (length(idxA) == 0 | length(idxB) == 0 | length(idxI) == 0)
    stop("Coefficient not found: ", varA, "/", varB, "/", inter_name)
  meanv <- c(b[idxA], b[idxB], b[idxI])
  covm  <- V[c(idxA, idxB, idxI), c(idxA, idxB, idxI)]
  reri    <- exp(sum(meanv)) - exp(meanv[1]) - exp(meanv[2]) + 1
  se_reri <- msm::deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                              mean = meanv, cov = covm)
  list(RERI = reri, se = se_reri,
       CI = c(reri - 1.96 * se_reri, reri + 1.96 * se_reri))
}
bootstrap_RERI <- function(data, formula, varA, varB,
                           inter_name = paste0(varA, ":", varB),
                           nboot = 1000, seed = 12345) {
  set.seed(seed); N <- nrow(data)
  out <- rep(NA_real_, nboot)
  for (i in seq_len(nboot)) {
    dboot <- data[sample(seq_len(N), replace = TRUE), , drop = FALSE]
    m <- tryCatch(coxph(formula, data = dboot), error = function(e) NULL)
    if (is.null(m)) next
    b <- coef(m)
    if (!all(c(varA, varB, inter_name) %in% names(b))) next
    out[i] <- exp(b[varA] + b[varB] + b[inter_name]) -
      exp(b[varA]) - exp(b[varB]) + 1
  }
  out <- out[!is.na(out)]
  list(est = mean(out), CI = quantile(out, c(.025, .975), na.rm = TRUE),
       n = length(out))
}
reri1 <- calc_RERI_from_cox(model_prs_high_afb, "PRS_high", "AFB_late")
reri1b <- bootstrap_RERI(dt_par, formula(model_prs_high_afb), "PRS_high", "AFB_late")
reri2 <- calc_RERI_from_cox(model_prs_afb, "PRS", "AFB_late")
reri2b <- bootstrap_RERI(dt_par, formula(model_prs_afb), "PRS", "AFB_late")
reri3 <- calc_RERI_from_cox(model_fh_par, "FH", "Parity_null")
reri3b <- bootstrap_RERI(dt, formula(model_fh_par), "FH", "Parity_null")
print(reri1); print(reri1b)
print(reri2); print(reri2b)
print(reri3); print(reri3b)
# 

# =====================================================================
# 5) main-effects（base / full）— 供 NRI/IDI/DCA/calibration/絕對風險
# =====================================================================
model_base <- coxph(
  Surv(baseline_age, end_age, status) ~ AFB_late + menarche_age_pref + FH +
    OC_ever_f + HRT_ever_f + menopause_status + oophorectomy_f +
    benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = analysis_afb, x = TRUE)
model_full <- coxph(
  Surv(baseline_age, end_age, status) ~ PRS + AFB_late + menarche_age_pref + FH +
    OC_ever_f + HRT_ever_f + menopause_status + oophorectomy_f +
    benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = analysis_afb, x = TRUE)
model_parity_main <- coxph(
  Surv(baseline_age, end_age, status) ~ PRS + Parity_null + menarche_age_pref + FH +
    OC_ever_f + HRT_ever_f + menopause_status + oophorectomy_f +
    benign_f + hysterectomy_ever_f +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = analysis_par)
C_base <- summary(model_base)$concordance[1]
C_full <- summary(model_full)$concordance[1]
cat("C base =", round(C_base,3), "| C full =", round(C_full,3),
    "| dC =", round(C_full - C_base,3), "\n")
# =====================================================================
# 6) 10-year risk / NRI / IDI / bootstrap
# =====================================================================
make_H0 <- function(model) {
  bh <- basehaz(model, centered = TRUE)
  t0 <- c(0, bh$time); h0 <- c(0, bh$hazard)
  function(t) approx(t0, h0, xout = t, method = "constant", f = 0, rule = 2)$y
}
risk10 <- function(model, newdata) {
  bh <- basehaz(model, centered = TRUE)
  lp <- predict(model, newdata = newdata, type = "lp")
  t0 <- c(0, bh$time); h0 <- c(0, bh$hazard)
  H0 <- function(t) approx(t0, h0, xout = t, method = "constant", f = 0, rule = 2)$y
  age0 <- newdata$baseline_age
  dH <- H0(age0 + 10) - H0(age0)
  1 - exp(-dH * exp(lp))
}
analysis_afb$risk_base <- risk10(model_base, analysis_afb)
analysis_afb$risk_prs  <- risk10(model_full, analysis_afb)
analysis_afb <- analysis_afb[!is.na(analysis_afb$risk_base) & !is.na(analysis_afb$risk_prs), ]
analysis_afb$follow <- analysis_afb$end_age - analysis_afb$baseline_age
analysis_afb$status_10y <- ifelse(analysis_afb$status == 1L & analysis_afb$follow <= 10, 1L, 0L)
keep10 <- (analysis_afb$status == 1L & analysis_afb$follow <= 10) | (analysis_afb$follow >= 10)
analysis_10y <- analysis_afb[keep10, ]
cat("10-year complete-case N =", nrow(analysis_10y),
    "| events =", sum(analysis_10y$status_10y), "\n")
cat_nri_pt <- function(p0, p1, y, cutoff) {
  c0 <- ifelse(p0 >= cutoff, 1L, 0L); c1 <- ifelse(p1 >= cutoff, 1L, 0L)
  up <- c1 > c0; down <- c1 < c0
  ev <- y == 1L; nev <- y == 0L
  (sum(up & ev) - sum(down & ev))/sum(ev) +
    (sum(down & nev) - sum(up & nev))/sum(nev)
}
cont_nri_pt <- function(p0, p1, y) {
  d <- p1 - p0; ev <- y == 1L; nev <- y == 0L
  (sum(d > 0 & ev)/sum(ev) - sum(d < 0 & ev)/sum(ev)) +
    (sum(d < 0 & nev)/sum(nev) - sum(d > 0 & nev)/sum(nev))
}
idi_pt <- function(p0, p1, y) {
  ev <- y == 1L; nev <- y == 0L
  (mean(p1[ev]) - mean(p0[ev])) - (mean(p1[nev]) - mean(p0[nev]))
}
boot_metrics <- function(data, nboot = 1000, seed = 123) {
  set.seed(seed); N <- nrow(data)
  res <- matrix(NA_real_, nboot, 4,
                dimnames = list(NULL, c("NRI_2%","NRI_3%","cNRI","IDI")))
  for (i in seq_len(nboot)) {
    dboot <- data[sample(seq_len(N), N, replace = TRUE), ]
    m0 <- tryCatch(coxph(Surv(baseline_age, end_age, status) ~
                           AFB_late + menarche_age_pref + FH +
                           OC_ever_f + HRT_ever_f + menopause_status + oophorectomy_f +
                           benign_f + hysterectomy_ever_f +
                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                         data = dboot), error = function(e) NULL)
    m1 <- tryCatch(coxph(Surv(baseline_age, end_age, status) ~
                           PRS + AFB_late + menarche_age_pref + FH +
                           OC_ever_f + HRT_ever_f + menopause_status + oophorectomy_f +
                           benign_f + hysterectomy_ever_f +
                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                         data = dboot), error = function(e) NULL)
    if (is.null(m0) || is.null(m1)) next
    r0 <- risk10(m0, dboot); r1 <- risk10(m1, dboot)
    ok <- !is.na(r0) & !is.na(r1)
    if (sum(ok) < 20) next
    follow <- dboot$end_age - dboot$baseline_age
    status10 <- ifelse(dboot$status == 1L & follow <= 10, 1L, 0L)
    keep <- (dboot$status == 1L & follow <= 10) | (follow >= 10)
    y <- status10[keep]; p0 <- r0[keep]; p1 <- r1[keep]
    if (sum(y) < 1 || sum(1 - y) < 1) next
    res[i, "NRI_2%"] <- cat_nri_pt(p0, p1, y, 0.02)
    res[i, "NRI_3%"] <- cat_nri_pt(p0, p1, y, 0.03)
    res[i, "cNRI"]   <- cont_nri_pt(p0, p1, y)
    res[i, "IDI"]    <- idi_pt(p0, p1, y)
  }
  lapply(colnames(res), function(k) {
    v <- na.omit(res[, k])
    c(est = mean(v), lo = quantile(v, .025, na.rm = TRUE),
      hi = quantile(v, .975, na.rm = TRUE), n_valid = length(v))
  })
}
bs <- boot_metrics(analysis_afb, nboot = 1000, seed = 123)
print(bs)
# 
# =====================================================================
# 7) DCA + ΔNB（paired bootstrap）
# =====================================================================
nb_at <- function(p, y, pt) {
  pred <- p >= pt
  tp <- sum(pred & y == 1L); fp <- sum(pred & y == 0L); n <- length(y)
  (tp - fp * pt/(1 - pt)) / n
}
thresholds <- seq(0, 0.10, by = 0.01)
nb_base <- sapply(thresholds, function(pt) nb_at(analysis_10y$risk_base, analysis_10y$status_10y, pt))
nb_prs  <- sapply(thresholds, function(pt) nb_at(analysis_10y$risk_prs,  analysis_10y$status_10y, pt))
prev <- mean(analysis_10y$status_10y)
nb_all <- prev - (1 - prev) * thresholds/(1 - thresholds)
deltaNB_2 <- nb_prs[which.min(abs(thresholds - 0.02))] - nb_base[which.min(abs(thresholds - 0.02))]
deltaNB_3 <- nb_prs[which.min(abs(thresholds - 0.03))] - nb_base[which.min(abs(thresholds - 0.03))]
nboot <- 1000; set.seed(123)
deltaNB_boot <- replicate(nboot, {
  db <- analysis_10y[sample(nrow(analysis_10y), replace = TRUE), ]
  nb_at(db$risk_prs, db$status_10y, 0.02) - nb_at(db$risk_base, db$status_10y, 0.02)
})
deltaNB3_boot <- replicate(nboot, {
  db <- analysis_10y[sample(nrow(analysis_10y), replace = TRUE), ]
  nb_at(db$risk_prs, db$status_10y, 0.03) - nb_at(db$risk_base, db$status_10y, 0.03)
})
cat("delta NB at 2% =", round(deltaNB_2,5),
    "| 95% CI:", quantile(deltaNB_boot, c(.025,.975), na.rm = TRUE), "\n")
cat("delta NB at 3% =", round(deltaNB_3,5),
    "| 95% CI:", quantile(deltaNB3_boot, c(.025,.975), na.rm = TRUE), "\n")
png("DCA_plot.png", width = 2000, height = 1600, res = 300)
plot(thresholds, nb_prs, type = "l", col = "blue", lwd = 2,
     xlab = "Threshold probability", ylab = "Net benefit",
     ylim = range(c(nb_base, nb_prs, nb_all, 0), na.rm = TRUE))
lines(thresholds, nb_base, col = "red", lwd = 2)
lines(thresholds, nb_all, col = "grey50", lty = 2)
abline(h = 0, lty = 3)
legend("topright",
       legend = c("PRS-integrated model", "Clinical-reproductive-familial base model", "Treat all", "Treat none"),
       col = c("blue","red","grey50","black"), lty = c(1,1,2,3), lwd = 2)
dev.off()
# 
#=====================================================================
  # 8) Calibration（logistic on 10-year outcome）
  # =====================================================================
analysis_10y$logit_risk <- qlogis(pmin(pmax(analysis_10y$risk_prs, 1e-6), 1 - 1e-6))
mean_pred <- mean(analysis_10y$risk_prs)
mean_obs  <- mean(analysis_10y$status_10y)
cat("mean predicted 10y risk =", round(mean_pred, 4),
    "\nobserved 10y risk      =", round(mean_obs, 4),
    "\ncalibration-in-the-large =", round(mean_pred - mean_obs, 4), "\n")
cal <- glm(status_10y ~ logit_risk, data = analysis_10y, family = binomial)
summary(cal)
nboot <- 1000; set.seed(123)
slopes <- replicate(nboot, {
  db <- analysis_10y[sample(nrow(analysis_10y), replace = TRUE), ]
  m <- glm(status_10y ~ logit_risk, data = db, family = binomial)
  coef(m)["logit_risk"]
})
cat("Calibration slope:", coef(cal)["logit_risk"],
    "\n95% CI:", quantile(slopes, c(.025, .975), na.rm = TRUE), "\n")
analysis_10y$decile <- ntile(analysis_10y$risk_prs, 10)
calibration_table <- analysis_10y %>%
  group_by(decile) %>%
  summarise(Predicted = mean(risk_prs), Observed = mean(status_10y), n = n())
print(calibration_table)
ggplot(calibration_table, aes(Predicted, Observed)) +
  geom_point(size = 3, colour = "blue") +
  geom_line(linewidth = 0.8, colour = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
  xlim(0, max(calibration_table$Predicted)) +
  ylim(0, max(calibration_table$Observed)) +
  labs(x = "Mean Predicted 10-year Risk", y = "Observed 10-year Risk") +
  theme_classic(base_size = 14)
ggsave("Calibration_plot.png", width = 6, height = 5, dpi = 600)
# 
#=====================================================================
  # 9) 絕對風險主分析 = reference profile（Figure 2 / crossing / 固定年齡）
  # =====================================================================
make_ref <- function(df) {
  ref <- as.data.frame(df[1, , drop = FALSE])
  for (nm in names(ref)) {
    if (is.factor(ref[[nm]]))
      ref[[nm]] <- factor(levels(ref[[nm]])[1], levels = levels(ref[[nm]]))
    else if (is.numeric(ref[[nm]])) ref[[nm]] <- 0
  }
  ref$menarche_age_pref <- 13
  ref
}
ref <- make_ref(analysis_afb)
ref_par <- make_ref(analysis_par)
H0      <- make_H0(model_full)
H0_par  <- make_H0(model_parity_main)
risk_ref <- function(model, H0fun, ref, prs_value, ages) {
  ref$PRS <- prs_value
  lp <- predict(model, newdata = ref, type = "lp")
  sapply(ages, function(a) 1 - exp(-(H0fun(a + 10) - H0fun(a)) * exp(lp)))
}
q <- quantile(analysis_afb$PRS, probs = c(.20, .80, .90, .95), na.rm = TRUE)
prs_values <- c(
  "Reference 20-80%" = mean(analysis_afb$PRS[analysis_afb$PRS >= q[1] & analysis_afb$PRS < q[2]], na.rm = TRUE),
  "Top 20%"          = mean(analysis_afb$PRS[analysis_afb$PRS >= q[2]], na.rm = TRUE),
  "Top 10%"          = mean(analysis_afb$PRS[analysis_afb$PRS >= q[3]], na.rm = TRUE),
  "Top 5%"           = mean(analysis_afb$PRS[analysis_afb$PRS >= q[4]], na.rm = TRUE)
)
ages <- seq(35, 70, by = 0.5)     # 資料支援範圍；<40 為外推
# 固定年齡表（40/50/60，主分析數值）
fix_ages <- c(40, 50, 60)
fix_tab <- do.call(rbind, lapply(names(prs_values), function(g) {
  data.frame(group = g,
             t(setNames(risk_ref(model_full, H0, ref, prs_values[[g]], fix_ages),
                        paste0("age_", fix_ages))))
}))
print(fix_tab)
# crossing（2% 與 3%，同一個 risk_ref）
get_age_cross <- function(age, risk, threshold) {
  o <- order(age); age <- age[o]; risk <- risk[o]
  idx <- which(risk[-1] >= threshold & risk[-length(risk)] < threshold)[1]
  if (is.na(idx)) return(NA_real_)
  age[idx] + (threshold - risk[idx]) / (risk[idx+1] - risk[idx]) * (age[idx+1] - age[idx])
}
age_cross_2 <- sapply(names(prs_values), function(g)
  get_age_cross(ages, risk_ref(model_full, H0, ref, prs_values[[g]], ages), 0.02))
age_cross_3 <- sapply(names(prs_values), function(g)
  get_age_cross(ages, risk_ref(model_full, H0, ref, prs_values[[g]], ages), 0.03))
print(age_cross_2); print(age_cross_3)
# bootstrap CI（定稿 nboot = 1000）
nboot <- 100; set.seed(123); N <- nrow(analysis_afb)
boot_risk <- array(NA_real_, dim = c(nboot, length(prs_values), length(ages)),
                   dimnames = list(NULL, names(prs_values), NULL))
for (b in seq_len(nboot)) {
  db <- analysis_afb[sample(seq_len(N), N, replace = TRUE), ]
  mb <- update(model_full, data = db)        # 與 model_full 同 formula
  H0b <- make_H0(mb)
  for (k in seq_along(prs_values))
    boot_risk[b, k, ] <- risk_ref(mb, H0b, ref, prs_values[[k]], ages)
}
ci_lo <- apply(boot_risk, c(2,3), function(x) quantile(x, .025, na.rm = TRUE))
ci_hi <- apply(boot_risk, c(2,3), function(x) quantile(x, .975, na.rm = TRUE))
cross_boot_2 <- apply(boot_risk, c(1,2), function(v) get_age_cross(ages, v, 0.02))
cross_boot_3 <- apply(boot_risk, c(1,2), function(v) get_age_cross(ages, v, 0.03))
print("2% crossing bootstrap 95% CI:")
print(apply(cross_boot_2, 2, function(x) quantile(x, c(.025,.975), na.rm = TRUE)))
print("3% crossing bootstrap 95% CI:")
print(apply(cross_boot_3, 2, function(x) quantile(x, c(.025,.975), na.rm = TRUE)))
risk_df <- do.call(rbind, lapply(names(prs_values), function(g)
  data.frame(group = g, age = ages,
             risk_10y = risk_ref(model_full, H0, ref, prs_values[[g]], ages))))
ci_df <- do.call(rbind, lapply(seq_along(prs_values), function(k)
  data.frame(group = names(prs_values)[k], age = ages, lo = ci_lo[k, ], hi = ci_hi[k, ])))
plot_df <- merge(risk_df, ci_df, by = c("group","age"))
plot_df$group <- factor(plot_df$group, levels = names(prs_values))
p_fig2 <- ggplot(plot_df, aes(age, risk_10y, color = group, fill = group)) +
  annotate("rect", xmin = 38, xmax = 40, ymin = -Inf, ymax = Inf,
           fill = "grey60", alpha = 0.25) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(0.02, 0.03), linetype = 2, color = "black") +
  annotate("text", x = 64.4, y = c(0.021, 0.031), label = c("2%","3%"),
           hjust = 0, size = 3.5) +
  annotate("text", x = 39, y = max(plot_df$hi, na.rm = TRUE) * 0.9,
           label = "Extrapolated\n(<40 y)", size = 3.2, color = "black") +
  scale_x_continuous(limits = c(38, 65.5), breaks = seq(40, 65, 5), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = "Age (years)", y = "10-year absolute breast cancer risk",
       color = "PRS group", fill = "PRS group") +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = NA))
ggsave("Figure2_reference_profile.png", p_fig2, width = 8, height = 5.5, dpi = 300)
# 
#=====================================================================
  # 10) 補充：marginal（Supplementary Table 6；35–70 歲每 5 歲）
  # =====================================================================
marginal_risk <- function(model, data, prs_value, ages, H0fun) {
  d <- as.data.frame(data); d$PRS <- prs_value
  lp <- predict(model, newdata = d, type = "lp")
  sapply(ages, function(a) {
    Sa   <- mean(exp(-H0fun(a)      * exp(lp)))
    Sa10 <- mean(exp(-H0fun(a + 10) * exp(lp)))
    1 - Sa10 / Sa
  })
}
ages_m <- seq(35, 70, by = 5)
marg_tab <- sapply(names(prs_values), function(g)
  marginal_risk(model_full, analysis_afb, prs_values[[g]], ages_m, H0))
rownames(marg_tab) <- ages_m
print(round(marg_tab * 100, 1))    # 填 Supplementary Table 6
# 
  
# =====================================================================
# 11) Figure 3：各因子曲線（reference profile，與主分析一致）
# =====================================================================
plot_cov_ref <- function(model, H0fun, ref, df, cov, values, ages, color_lab, file) {
  risk_df <- do.call(rbind, lapply(names(values), function(lab) {
    nd <- ref
    if (is.factor(df[[cov]])) {
      nd[[cov]] <- factor(values[[lab]], levels = levels(df[[cov]]))
    } else {
      nd[[cov]] <- as.numeric(values[[lab]])
    }
    lp <- predict(model, newdata = nd, type = "lp")
    data.frame(group = lab, age = ages,
               risk_10y = sapply(ages, function(a)
                 1 - exp(-(H0fun(a + 10) - H0fun(a)) * exp(lp))))
  }))
  p <- ggplot(risk_df, aes(age, risk_10y, color = group)) +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0.02, linetype = "dashed", color = "red") +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(x = "Age (years)", y = "10-year absolute risk", color = color_lab) +
    theme_classic(base_size = 13)
  ggsave(file, p, width = 5.2, height = 3.8, dpi = 300)
}
plot_cov_ref(model_full, H0, ref, analysis_afb, "FH",
             c("No family history" = 0, "With family history" = 1),
             ages, "Family history", "fig3_A_fh.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "AFB_late",
             c("First birth <30" = 0, "First birth \u226530" = 1),
             ages, "Age at first birth", "fig3_B_afb.png")
plot_cov_ref(model_parity_main, H0_par, ref_par, analysis_par, "Parity_null",
             c("Parous" = 0, "No live parity" = 1),
             ages, "Parity", "fig3_C_parity.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "menarche_age_pref",
             c("Early (11 y)" = 11, "Average (13 y)" = 13, "Late (15 y)" = 15),
             ages, "Age at menarche", "fig3_D_menarche.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "OC_ever_f",
             c("Never OC" = "No", "Ever OC" = "Yes"),
             ages, "Oral contraceptive use", "fig3_E_oc.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "HRT_ever_f",
             c("Never HRT" = "No", "Ever HRT" = "Yes"),
             ages, "HRT use", "fig3_F_hrt.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "menopause_status",
             c("Premenopausal" = "Premenopausal",
               "Postmenopausal" = "Postmenopausal",
               "Uncertain" = "Uncertain"),
             ages, "Menopausal status", "fig3_G_menopause.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "oophorectomy_f",
             c("No oophorectomy" = "No", "Bilateral oophorectomy" = "Yes"),
             ages, "Oophorectomy", "fig3_H_oophorectomy.png")
plot_cov_ref(model_full, H0, ref, analysis_afb, "benign_f",
             c("No benign disease" = "No", "Benign breast disease" = "Yes"),
             ages, "Benign breast disease", "fig3_I_benign.png")
# =====================================================================
# 12) Table 1（gtsummary；hysterectomy 改三級以配合模型）
# =====================================================================
library(gtsummary)
need <- c("PRS","baseline_age","menarche_age_pref","follow_up_years",
          "Parity_null","AFB_late","FH","status",
          "2784-0.0","X2814.0.0","2834-0.0","3591-0.0",
          "benign_breast_any_prebaseline","2724-0.0",
          "21001-0.0","20116-0.0","1558-0.0")
stopifnot(all(need %in% names(dt)))
dt_tbl <- dt %>%
  mutate(
    PRS_stratum = factor(case_when(
      PRS < quantile(PRS, 0.20, na.rm = TRUE) ~ "Lowest 20%",
      PRS < quantile(PRS, 0.80, na.rm = TRUE) ~ "Middle 20\u201380%",
      PRS < quantile(PRS, 0.90, na.rm = TRUE) ~ "80\u2013<90%",
      PRS < quantile(PRS, 0.95, na.rm = TRUE) ~ "90\u2013<95%",
      TRUE ~ "Highest 5%"),
      levels = c("Lowest 20%","Middle 20\u201380%","80\u2013<90%","90\u2013<95%","Highest 5%")),
    OC_ever = factor(case_when(`2784-0.0` == 0 ~ "No", `2784-0.0` == 1 ~ "Yes",
                               TRUE ~ NA_character_), levels = c("No","Yes")),
    HRT_ever = factor(case_when(X2814.0.0 == 0 ~ "No", X2814.0.0 == 1 ~ "Yes",
                                TRUE ~ NA_character_), levels = c("No","Yes")),
    oophorectomy = factor(case_when(`2834-0.0` == 0 ~ "No", `2834-0.0` == 1 ~ "Yes",
                                    TRUE ~ NA_character_), levels = c("No","Yes")),
    hysterectomy = factor(case_when(`3591-0.0` == 0 ~ "No",
                                    `3591-0.0` == 1 ~ "Yes",
                                    TRUE ~ "Uncertain"),
                          levels = c("No","Yes","Uncertain")),
    benign = factor(if_else(benign_breast_any_prebaseline, "Yes", "No"),
                    levels = c("No","Yes")),
    menopause = factor(case_when(`2724-0.0` == 0 ~ "Premenopausal",
                                 `2724-0.0` == 1 ~ "Postmenopausal",
                                 `2724-0.0` %in% c(2,3) ~ "Uncertain",
                                 TRUE ~ NA_character_),
                       levels = c("Premenopausal","Postmenopausal","Uncertain")),
    smoking = factor(case_when(`20116-0.0` == 0 ~ "Never",
                               `20116-0.0` == 1 ~ "Former",
                               `20116-0.0` == 2 ~ "Current",
                               TRUE ~ NA_character_),
                     levels = c("Never","Former","Current")),
    alcohol_ever = factor(case_when(`1558-0.0` %in% 1:5 ~ "Yes",
                                    `1558-0.0` == 6 ~ "No",
                                    TRUE ~ NA_character_), levels = c("No","Yes"))
  )
tbl1 <- dt_tbl %>%
  select(PRS_stratum, age = baseline_age, menarche = menarche_age_pref,
         BMI = `21001-0.0`, follow_up = follow_up_years,
         Parity_null, AFB_late, OC_ever, HRT_ever, menopause,
         oophorectomy, hysterectomy, benign, FH,
         smoking, alcohol_ever, status) %>%
  tbl_summary(
    by = PRS_stratum,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    missing = "no",
    digits = list(all_continuous() ~ 1, all_categorical() ~ 1)
  ) %>%
  add_overall() %>%
  add_p(test = list(all_continuous() ~ "kruskal.test",
                    all_categorical() ~ "chisq.test")) %>%
  modify_header(label = "**Characteristic**") %>%
  bold_labels()
tbl1
library(flextable)
tbl1 %>% as_flex_table() %>% save_as_docx(path = "Table1.docx")
# person-year / crude incidence（Table 1 末兩列）
incidence_summary <- dt_tbl %>%
  group_by(PRS_stratum) %>%
  summarise(n = n(), events = sum(status == 1, na.rm = TRUE),
            person_years = sum(follow_up_years, na.rm = TRUE),
            crude_incidence_rate = events / person_years * 1000,
            .groups = "drop")
print(incidence_summary)