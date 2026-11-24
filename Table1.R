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