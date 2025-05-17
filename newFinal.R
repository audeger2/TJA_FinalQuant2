library(tidyverse)
library(sensemakr)
library(rdrobust)
library(optmatch)
library(RItools)
library(janitor)
library(pwr)
library(fixest)
library(readxl)
library(gt)
library(tibble)
library(ggplot2)

#Load the dataset
df <- read_csv("/Users/andrewdeger/Downloads/final_data.csv")

#Create derived variables and define treatment
df <- df %>%
  mutate(
    abs_pvi_change = abs(`Change in PVI`),
    abs_change_cfdyn = abs(`Change in CFDyn`)
  ) %>%
  group_by(Congress) %>%
  mutate(
    quantile_cutoff = quantile(abs_pvi_change, 0.75, na.rm = TRUE),
    treat = ifelse(abs_pvi_change >= quantile_cutoff, 1, 0)
  ) %>%
  ungroup()

#Matching for 113th Congress
data_113 <- df %>%
  filter(Congress == 113, !is.na(treat), !is.na(abs_change_cfdyn),
         !is.na(`Vote Share`), !is.na(`Income (Median)`), !is.na(`CF score`))

distance_113 <- match_on(
  treat ~ `Vote Share` + `Income (Median)` + `CF score`,
  data = data_113,
  method = "mahalanobis"
)

pairs_113 <- pairmatch(distance_113, data = data_113)

matched_113 <- data_113 %>%
  mutate(pair_id = pairs_113) %>%
  filter(!is.na(pair_id))

bal_113 <- balanceTest(
  treat ~ `Vote Share` + `Income (Median)` + `CF score` + strata(pair_id),
  data = matched_113
)


#Calculate number of units not matched
#113th Congress
n_total_113 <- nrow(data_113)
n_matched_113 <- nrow(matched_113)
n_unmatched_113 <- n_total_113 - n_matched_113

cat("113th Congress:\n")
cat("  Total units before matching:", n_total_113, "\n")
cat("  Matched units:", n_matched_113, "\n")
cat("  Unmatched units:", n_unmatched_113, "\n\n")

#118th Congress
n_total_118 <- nrow(data_118)
n_matched_118 <- nrow(matched_118)
n_unmatched_118 <- n_total_118 - n_matched_118

cat("118th Congress:\n")
cat("  Total units before matching:", n_total_118, "\n")
cat("  Matched units:", n_matched_118, "\n")
cat("  Unmatched units:", n_unmatched_118, "\n")


#Format balance table
bal_df <- as.data.frame(bal_113$results) %>%
  rownames_to_column(var = "Covariate") %>%
  select(Covariate, adj.diff = `adj.diff.--`, std.diff = `std.diff.--`, z = `z.--`) %>%
  rename(`Adjusted Difference` = adj.diff, `Standardized Difference` = std.diff, `Z-Score` = z)

p_bal_113 <- gt(bal_df) %>%
  tab_header(title = "Table 1: Covariate Balance After Matching (113th Congress)")
print(p_bal_113)

model_att_113 <- lm(abs_change_cfdyn ~ treat, data = matched_113)
summary(model_att_113)

model_att_cov_113 <- lm(abs_change_cfdyn ~ treat + `CF score` + `Vote Share` + `Income (Median)`, data = matched_113)

sensitivity_113 <- sensemakr(
  model = model_att_cov_113,
  treatment = "treat",
  benchmark_covariates = "`CF score`",
  kd = 1:3
)
summary(sensitivity_113)
plot(sensitivity_113)

#Matching for 118th Congress
data_118 <- df %>%
  filter(Congress == 118, !is.na(treat), !is.na(abs_change_cfdyn),
         !is.na(`Vote Share`), !is.na(`Income (Median)`), !is.na(`CF score`))

distance_118 <- match_on(
  treat ~ `Vote Share` + `Income (Median)` + `CF score`,
  data = data_118,
  method = "mahalanobis"
)

pairs_118 <- pairmatch(distance_118, data = data_118)

matched_118 <- data_118 %>%
  mutate(pair_id = pairs_118) %>%
  filter(!is.na(pair_id))

bal_118 <- balanceTest(
  treat ~ `Vote Share` + `Income (Median)` + `CF score` + strata(pair_id),
  data = matched_118
)

bal_df_118 <- as.data.frame(bal_118$results) %>%
  rownames_to_column(var = "Covariate") %>%
  select(Covariate, adj.diff = `adj.diff.--`, std.diff = `std.diff.--`, z = `z.--`) %>%
  rename(`Adjusted Difference` = adj.diff, `Standardized Difference` = std.diff, `Z-Score` = z)

p_bal_118 <- gt(bal_df_118) %>%
  tab_header(title = "Table 2: Covariate Balance After Matching (118th Congress)")
print(p_bal_118)

model_att_118 <- lm(abs_change_cfdyn ~ treat, data = matched_118)
summary(model_att_118)

model_att_cov_118 <- lm(abs_change_cfdyn ~ treat + `CF score` + `Vote Share` + `Income (Median)`, data = matched_118)

sensitivity_118 <- sensemakr(
  model = model_att_cov_118,
  treatment = "treat",
  benchmark_covariates = "`CF score`",
  kd = 1:3
)
summary(sensitivity_118)
plot(sensitivity_118)

#Fixed Effects Regression
reg_data <- df %>%
  filter(
    !is.na(abs_change_cfdyn),
    !is.na(abs_pvi_change),
    !is.na(`Vote Share`),
    !is.na(`Income (Median)`),
    !is.na(`CF score`),
    !is.na(`Congress`)
  )

model_fe <- feols(
  abs_change_cfdyn ~ abs_pvi_change + `Vote Share` + `Income (Median)` + `CF score` | `Congress`,
  data = reg_data
)
summary(model_fe)

#Tidy the model
fe_tidy <- broom::tidy(model_fe)

#Reformat
fe_table <- fe_tidy %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  rename(
    Variable = term,
    `Coefficient` = estimate,
    `Std. Error` = std.error,
    `t-stat` = statistic,
    `p-value` = p.value
  ) %>%
  select(Variable, Coefficient, `Std. Error`, `t-stat`, `p-value`, Significance)

#Create the table
p_fe_table <- gt(fe_table) %>%
  tab_header(
    title = "Table 3: Fixed Effects Regression Results"
  ) %>%
  fmt_number(
    columns = c(Coefficient, `Std. Error`, `t-stat`, `p-value`),
    decimals = 3
  ) %>%
  cols_align(align = "center")

#Print the table
print(p_fe_table)


#Fixed Effects Scatterplot
#113th Congress
demeaned_113 <- reg_data %>%
  filter(Congress == 113) %>%
  mutate(
    y_within = abs_change_cfdyn - mean(abs_change_cfdyn, na.rm = TRUE),
    x_within = abs_pvi_change - mean(abs_pvi_change, na.rm = TRUE)
  )

p_fe_113 <- ggplot(demeaned_113, aes(x = x_within, y = y_within)) +
  geom_point(alpha = 0.5, color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Within-Congress Abs Change in PVI",
    y = "Within-Congress Abs Change in CF Score"
  ) +
  theme_minimal(base_size = 13)

print(p_fe_113)

#118th Congress
demeaned_118 <- reg_data %>%
  filter(Congress == 118) %>%
  mutate(
    y_within = abs_change_cfdyn - mean(abs_change_cfdyn, na.rm = TRUE),
    x_within = abs_pvi_change - mean(abs_pvi_change, na.rm = TRUE)
  )

p_fe_118 <- ggplot(demeaned_118, aes(x = x_within, y = y_within)) +
  geom_point(alpha = 0.5, color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Within-Congress Abs Change in PVI",
    y = "Within-Congress Abs Change in CF Score"
  ) +
  theme_minimal(base_size = 13)

print(p_fe_118)


#Powe analysis/test
sd_113 <- sd(matched_113$abs_change_cfdyn, na.rm = TRUE)
n_113 <- nrow(matched_113) / 2
mde_113 <- pwr.t.test(n = n_113, power = 0.80, sig.level = 0.05, type = "two.sample")
mde_raw_113 <- mde_113$d * sd_113

cat("113th Congress:\n")
cat("  SD of outcome:", round(sd_113, 4), "\n")
cat("  MDE (Cohen's d):", round(mde_113$d, 4), "\n")
cat("  MDE in raw units:", round(mde_raw_113, 4), "\n\n")

sd_118 <- sd(matched_118$abs_change_cfdyn, na.rm = TRUE)
n_118 <- nrow(matched_118) / 2
mde_118 <- pwr.t.test(n = n_118, power = 0.80, sig.level = 0.05, type = "two.sample")
mde_raw_118 <- mde_118$d * sd_118

cat("118th Congress:\n")
cat("  SD of outcome:", round(sd_118, 4), "\n")
cat("  MDE (Cohen's d):", round(mde_118$d, 4), "\n")
cat("  MDE in raw units:", round(mde_raw_118, 4), "\n")


#Placebo tests
#113th Congress
set.seed(123)
n_sim <- 1000
n_treat <- sum(matched_113$treat == 1)
n_control <- sum(matched_113$treat == 0)
placebo_pvals <- replicate(n_sim, {
  placebo_treat <- sample(c(rep(1, n_treat), rep(0, n_control)))
  model <- lm(abs_change_cfdyn ~ placebo_treat + `Vote Share` + `Income (Median)` + `CF score`, data = matched_113)
  coef(summary(model))["placebo_treat", "Pr(>|t|)"]
})
false_positive_rate <- mean(placebo_pvals < 0.05)

p113 <- ggplot(tibble(p_value = placebo_pvals), aes(x = p_value)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    subtitle = paste0("Empirical false positive rate: ", round(false_positive_rate, 3)),
    x = "Placebo P-value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 13)
print(p113)

#118th Congress
set.seed(123)
n_treat_118 <- sum(matched_118$treat == 1)
n_control_118 <- sum(matched_118$treat == 0)
placebo_pvals_118 <- replicate(n_sim, {
  placebo_treat <- sample(c(rep(1, n_treat_118), rep(0, n_control_118)))
  model <- lm(abs_change_cfdyn ~ placebo_treat + `Vote Share` + `Income (Median)` + `CF score`, data = matched_118)
  coef(summary(model))["placebo_treat", "Pr(>|t|)"]
})
false_positive_rate_118 <- mean(placebo_pvals_118 < 0.05)

p118 <- ggplot(tibble(p_value = placebo_pvals_118), aes(x = p_value)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "white", boundary = 0) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    subtitle = paste0("Empirical false positive rate: ", round(false_positive_rate_118, 3)),
    x = "P-values from placebo treatment",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 13)
print(p118)


#Regression plots
#Filter and clean data for each Congress
regression_113 <- df %>%
  filter(Congress == 113, !is.na(abs_change_cfdyn), !is.na(abs_pvi_change))

regression_118 <- df %>%
  filter(Congress == 118, !is.na(abs_change_cfdyn), !is.na(abs_pvi_change))

#Plot for 113th Congress
p113_reg <- ggplot(regression_113, aes(x = abs_pvi_change, y = abs_change_cfdyn)) +
  geom_point(alpha = 0.6, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Absolute Change in District PVI",
    y = "Absolute Change in Legislator CF Score"
  ) +
  theme_minimal(base_size = 13)

#Print the 113th Congress plot
print(p113_reg)

#Plot for 118th Congress
p118_reg <- ggplot(regression_118, aes(x = abs_pvi_change, y = abs_change_cfdyn)) +
  geom_point(alpha = 0.6, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Absolute Change in District PVI",
    y = "Absolute Change in Legislator CF Score"
  ) +
  theme_minimal(base_size = 13)

# Print the 118th Congress plot
print(p118_reg)

