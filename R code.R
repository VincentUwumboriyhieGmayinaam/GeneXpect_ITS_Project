# Install packages required for the analysis (uncomment if needed)
install.packages("lmtest") ; install.packages("Epi")
install.packages("tsModel"); install.packages("vcd")
install.packages("readxl")
# load the packages
library(foreign) ; library(tsModel) ; library("lmtest") ; library("Epi")
library("splines") ; library("vcd"); library(readxl)

# read data from csv file
data <- read_excel("PAUL_MPH TB DATA_4ANALYSIS NEW (1).xlsx")
head(data)
View(data)

names(data) <- trimws(names(data))
names(data)
data$Month <- match(tolower(data$Month),
                    tolower(month.name))
data$time <- with(data, (Year - min(Year)) * 12 + Month)
summary(data$time)

library(dplyr)

monthly_tb <- data %>%
  group_by(Year, Month, time, GeneXpect) %>%
  summarise(
    tests = n(),
    positives = sum(Result == 1, na.rm = TRUE),
    detection_rate = (positives / tests) * 100,
    .groups = "drop")

head(monthly_tb)
summary(monthly_tb$detection_rate)

intervention_time <- min(monthly_tb$time[monthly_tb$GeneXpect == 1])
intervention_time

pre  <- monthly_tb[monthly_tb$time < intervention_time, ]
post <- monthly_tb[monthly_tb$time >= intervention_time, ]

fit_pre  <- lm(detection_rate ~ time, data = pre)
fit_post <- lm(detection_rate ~ time, data = post)


plot(monthly_tb$time,
     monthly_tb$detection_rate,
     type = "n",
     ylim = c(0, max(monthly_tb$detection_rate, na.rm = TRUE)),
     xlab = "Year",
     ylab = "TB detection rate (%)",
     bty = "l",
     xaxt = "n")

# rect(intervention_time,
#      0,
#      max(monthly_tb$time),
#      max(monthly_tb$detection_rate, na.rm = TRUE),
#      col = rgb(1.0, 1.0, 1.0, 1.0),  # light grey with transparency
#      border =grey(0.5, 0.5),
#      lwd = 1)

abline(v = intervention_time,
       col = "black",
       lwd = 2,
       lty = 2)

# Observed points
points(monthly_tb$time,
       monthly_tb$detection_rate, type = "l",
       pch = 16,
       cex = 0.6,
       col = "steelblue")

# Pre- and post-intervention linear trends
abline(fit_pre,  lwd = 2, col = "navy")
abline(fit_post, lwd = 2, col = "firebrick")

# X-axis (calendar years)
axis(1,
     at = seq(6, max(monthly_tb$time), by = 12),
     labels = sort(unique(monthly_tb$Year)))

# Title
title("TB Case Detection Trends Before and After GeneXpert Introduction")


legend("topleft",
       legend = c("Observed TB detection rate",
                  "Pre-intervention trend",
                  "Post-intervention trend",
                  "Post-intervention period"),
       col = c("steelblue",
               "navy",
               "firebrick",
               "black",
       pch = c(16, NA, NA, 15),
       lty = c(NA, 2, 2, NA),
       lwd = c(NA, 2, 2, NA),
       pt.cex = c(0.8, NA, NA, 1),
       cex = 0.8,
       bty = "n")

text(x = intervention_time + 
       (max(monthly_tb$time) - intervention_time) / 1.7,
     y = max(monthly_tb$detection_rate, na.rm = TRUE) * 0.9,
     labels = "Post-GeneXpert introduction",
     col = "black",
     cex = 0.9,
     font = 2)

# It is also useful to produce summary statistics
summary(data)

#tabulate aces before and after the gene-expect introduction
summary(data$Result[data$GeneXpect==0])
summary(data$Result[data$GeneXpect==1])

summary(monthly_tb$detection_rate[monthly_tb$GeneXpect==0])
summary(monthly_tb$detection_rate[monthly_tb$GeneXpect==1])


tb_model <- monthly_tb %>%
  mutate(
    intervention = GeneXpect,
    month = Month
  )

#Poisson with the standardised population as an offset
model1 <- glm(positives ~ offset(log(tests)) + intervention + time,
    family = poisson,
    data = tb_model)

summary(model1)
summary(model1)$dispersion
round(ci.lin(model1,Exp=T),3)

model2<-glm(positives ~ offset(log(tests)) + intervention + time,
    family = quasipoisson,
    data = tb_model)

summary(model2)
summary(model2)$dispersion
round(ci.lin(model2,Exp=T),3)

model3<-glm(positives ~ offset(log(tests)) + intervention + time +
      harmonic(month, 2, 12),
    family = quasipoisson,
    data = tb_model)

summary(model3)
summary(model3)$dispersion
round(ci.lin(model3,Exp=T),3)


model_step <- glm(
  positives ~ offset(log(tests)) + intervention + time +
    harmonic(month, 2, 12),
  family = quasipoisson,
  data = tb_model
)
summary(model_step)

model_slope <- glm(
  positives ~ offset(log(tests)) + intervention * time +
    harmonic(month, 2, 12),
  family = quasipoisson,
  data = tb_model
)
summary(model_slope)


anova(model_step, model_slope, test = "F")

# Prediction dataset (real-world scenario)
pred_data <- tb_model


# Fitted predictions (with intervention)
pred_data$fitted <- predict(
  model_slope,
  type = "response"
) / pred_data$tests * 100

# Counterfactual dataset (no GeneXpert ever)
cf_data <- pred_data
cf_data$intervention <- 0

cf_data$counterfactual <- predict(
  model_slope,
  newdata = cf_data,
  type = "response"
) / cf_data$tests * 100

pred_data$observed <- pred_data$detection_rate

plot_df <- pred_data %>%
  select(time, Year, observed, fitted) %>%
  left_join(
    cf_data %>% select(time, counterfactual),
    by = "time"
  )

library(ggplot2)

ymax <- max(plot_df$observed, na.rm = TRUE)
xmax <- max(plot_df$time)

ggplot(plot_df, aes(x = time)) +
  
  # Post-intervention window
  annotate("rect",
           xmin = intervention_time,
           xmax = xmax,
           ymin = 0,
           ymax = ymax,
           fill = "grey80",
           alpha = 0.5) +
  
  # Observed data
  geom_point(aes(y = observed, colour = "Observed"),
             size = 1.6) +
  
  # Fitted values
  geom_line(aes(y = fitted, colour = "Fitted"),
            linewidth = 1.2) +
  
  # Counterfactual
  geom_line(aes(y = counterfactual, colour = "Counterfactual"),
            linewidth = 1.2,
            linetype = "dashed") +
  
  # Intervention line
  geom_vline(xintercept = intervention_time,
             linetype = "dashed",
             colour = "grey40") +
  
  # X-axis as years
  scale_x_continuous(
    breaks = seq(6, xmax, by = 12),
    labels = sort(unique(plot_df$Year))
  ) +
  
  scale_colour_manual(
    values = c(
      "Observed" = "steelblue",
      "Fitted" = "navy",
      "Counterfactual" = "firebrick"
    )
  ) +
  
  labs(
    title = "Interrupted Time Series of TB Detection Following GeneXpert Introduction",
    x = "Year",
    y = "TB detection rate (%)",
    colour = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

impact_df <- plot_df %>%
  filter(time >= intervention_time)

impact_df <- impact_df %>%
  mutate(
    abs_impact = fitted - counterfactual
  )


impact_df <- impact_df %>%
  mutate(
    rel_impact = (fitted - counterfactual) / counterfactual * 100
  )

impact_summary <- impact_df %>%
  summarise(
    mean_absolute_impact = mean(abs_impact, na.rm = TRUE),
    mean_relative_impact = mean(rel_impact, na.rm = TRUE),
    median_absolute_impact = median(abs_impact, na.rm = TRUE),
    median_relative_impact = median(rel_impact, na.rm = TRUE)
  )

impact_summary

impact_12m <- impact_df %>%
  filter(time == intervention_time + 12) %>%
  select(abs_impact, rel_impact)

impact_12m

impact_df <- impact_df %>%
  mutate(
    cumulative_abs_impact = cumsum(abs_impact)
  )


library(MASS)

set.seed(123)

# Number of simulations
nsim <- 1000

# Extract coefficients and covariance matrix
beta_hat <- coef(model_slope)
V_beta   <- vcov(model_slope)

# Simulate coefficients
beta_sim <- mvrnorm(nsim, beta_hat, V_beta)


# Design matrix (observed scenario)
X_obs <- model.matrix(model_slope, data = pred_data)

# Design matrix (counterfactual: intervention = 0)
cf_data <- pred_data
cf_data$intervention <- 0
X_cf <- model.matrix(model_slope, data = cf_data)

# Simulated fitted predictions
mu_obs <- exp(X_obs %*% t(beta_sim)) / pred_data$tests * 100

# Simulated counterfactual predictions
mu_cf  <- exp(X_cf  %*% t(beta_sim)) / pred_data$tests * 100

impact_sim <- mu_obs - mu_cf

post_idx <- pred_data$time >= intervention_time

impact_ci <- apply(impact_sim[post_idx, ], 1, function(x) {
  c(
    mean = mean(x),
    lci  = quantile(x, 0.025),
    uci  = quantile(x, 0.975)
  )
})

impact_ci <- as.data.frame(t(impact_ci))
impact_ci$time <- pred_data$time[post_idx]

rel_impact_sim <- (mu_obs - mu_cf) / mu_cf * 100

rel_impact_ci <- apply(rel_impact_sim[post_idx, ], 1, function(x) {
  c(
    mean = mean(x),
    lci  = quantile(x, 0.025),
    uci  = quantile(x, 0.975)
  )
})

rel_impact_ci <- as.data.frame(t(rel_impact_ci))
rel_impact_ci$time <- pred_data$time[post_idx]


rel_impact_ci <- as.data.frame(t(rel_impact_ci))
str(rel_impact_ci)
impact_ci <- as.data.frame(t(impact_ci))


str(impact_ci)
str(rel_impact_ci)
impact_ci_fixed <- as.data.frame(t(impact_ci))

# Name the columns correctly
colnames(impact_ci_fixed) <- c("mean", "lci", "uci", "time")

# Convert to numeric (important)
impact_ci_fixed[] <- lapply(impact_ci_fixed, as.numeric)

str(impact_ci_fixed)


impact_table <- data.frame(
  Measure = c(
    "Mean absolute impact (pp)",
    "Median absolute impact (pp)",
    "Mean relative impact (%)",
    "Median relative impact (%)"
  ),
  
  Estimate = c(
    mean(impact_ci_fixed$mean, na.rm = TRUE),
    median(impact_ci_fixed$mean, na.rm = TRUE),
    mean(rel_impact_ci$mean, na.rm = TRUE),
    median(rel_impact_ci$mean, na.rm = TRUE)
  ),
  
  Lower_95CI = c(
    mean(impact_ci_fixed$lci, na.rm = TRUE),
    median(impact_ci_fixed$lci, na.rm = TRUE),
    mean(rel_impact_ci$lci, na.rm = TRUE),
    median(rel_impact_ci$lci, na.rm = TRUE)
  ),
  
  Upper_95CI = c(
    mean(impact_ci_fixed$uci, na.rm = TRUE),
    median(impact_ci_fixed$uci, na.rm = TRUE),
    mean(rel_impact_ci$uci, na.rm = TRUE),
    median(rel_impact_ci$uci, na.rm = TRUE)
  )
)

round(impact_table, 2)

impact_table

res2 <- residuals(model2,type="deviance")
plot(data$time,res2,ylim=c(-5,10),pch=19,cex=0.7,col=grey(0.6),
     main="Residuals over time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2)

acf(res2)
pacf(res2)

model3 <- glm(aces ~ offset(log(stdpop)) + smokban + time + 
                harmonic(month,2,12), family=quasipoisson, data)
summary(model3)
summary(model3)$dispersion
round(ci.lin(model3,Exp=T),3)

# EFFECTS
ci.lin(model3,Exp=T)["smokban",5:7]

# TREND
exp(coef(model3)["time"]*12)


res3 <- residuals(model3,type="deviance")
plot(res3,ylim=c(-5,10),pch=19,cex=0.7,col=grey(0.6),main="Residuals over time",
     ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2)
acf(res3)
pacf(res3)

# not surprisingly, the p-value is similar to that of the interaction term

# 

