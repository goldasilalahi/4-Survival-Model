# Import libraries
library(survival)
library(KMsurv)
library(MASS)
library(survminer)
library(readxl)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(car)
library(gridExtra)

# Load data
data <- read_excel("/Users/Golda/Github/College/sem 4/Model Survival/Final Project/melanoma.xlsx")

# Remove rows with status 3
data <- data %>% filter(status != 3)

# Check data structure
head(data)
str(data)
dim(data)
summary(data)

# Convert categorical variables
data$sex <- as.factor(data$sex)
data$ulcer <- as.factor(data$ulcer)

# Convert numerical variables to categorical
data$thickness_strata <- cut(data$thickness, 
                             breaks = c(-Inf, 1, 2, 4, Inf), 
                             labels = c("Level 1", "Level 2", "Level 3", "Level 4"),
                             right = FALSE)

data$age_category <- cut(data$age, 
                         breaks = c(4, 11, 25, 45, 65, Inf), 
                         labels = c("kanak-kanak", "remaja", "dewasa", "lansia", "manula"),
                         right = FALSE)

# Define plot variables, types, and colors
plot_vars <- c("time", "status", "sex", "age_category", "year", "thickness_strata", "ulcer")
plot_types <- c("hist", "bar", "bar", "bar", "hist", "bar", "bar")
plot_colors <- c("steelblue", "lightcoral", "maroon", "mediumseagreen", "gold", "turquoise", "lightsalmon")

# Generate plots using a loop
plots <- lapply(seq_along(plot_vars), function(i) {
  var <- plot_vars[i]
  p <- ggplot(data, aes_string(x = var)) + ggtitle(paste("Distribution of", gsub("_", " ", var)))
  
  if (plot_types[i] == "hist") {
    p + geom_histogram(binwidth = ifelse(var == "year", 1, 50), fill = plot_colors[i], color = "black")
  } else {
    p + geom_bar(fill = plot_colors[i], color = "black")
  }
})

# Display all plots
library(gridExtra)
grid.arrange(grobs = plots, ncol = 2)

# Kaplan-Meier plots
St.ulcer <- surv_fit(Surv(time, status == 1) ~ ulcer, data = data)
ggsurvplot(St.ulcer, data = data, conf.int = TRUE, pval = TRUE, risk.table = TRUE, risk.table.height = .3)

St.sex <- surv_fit(Surv(time, status == 1) ~ sex, data = data)
ggsurvplot(St.sex, data = data, conf.int = TRUE, pval = TRUE, risk.table = TRUE, risk.table.height = .3)

St.age_category <- surv_fit(Surv(time, status == 1) ~ age_category, data = data)
ggsurvplot(St.age_category, data = data, conf.int = TRUE, pval = TRUE, risk.table = TRUE, risk.table.height = .3)

St.thickness_strata <- surv_fit(Surv(time, status == 1) ~ thickness_strata, data = data)
ggsurvplot(St.thickness_strata, data = data, conf.int = TRUE, pval = TRUE, risk.table = TRUE, risk.table.height = .3)

# Log-rank test
surv_object <- Surv(time = data$time, event = data$status == 1)

logrank_sex <- survdiff(surv_object ~ sex, data = data)
logrank_ulcer <- survdiff(surv_object ~ ulcer, data = data)
logrank_thickness_strata <- survdiff(surv_object ~ thickness_strata, data = data)
logrank_age_category <- survdiff(surv_object ~ age_category, data = data)

logrank_sex
logrank_ulcer
logrank_thickness_strata
logrank_age_category

# Cox proportional hazards model
model_ulcer <- coxph(Surv(time, status) ~ ulcer, data = data)
model_sex <- coxph(Surv(time, status) ~ sex, data = data)
model_age_category <- coxph(Surv(time, status) ~ age_category, data = data)
model_thickness_strata <- coxph(Surv(time, status) ~ thickness_strata, data = data)

model_ulcer
model_sex
model_age_category
model_thickness_strata

ggforest(model_ulcer, data = data)
ggforest(model_sex, data = data)
ggforest(model_age_category, data = data)
ggforest(model_thickness_strata, data = data)

# Model selection
full_model <- coxph(surv_object ~ sex + age_category + year + thickness_strata + ulcer, data = data)
summary(full_model)

both_model <- stepAIC(full_model, direction = "both", trace = TRUE)
summary(both_model)

forward_model <- stepAIC(full_model, direction = "forward", trace = TRUE)
summary(forward_model)

backward_model <- stepAIC(full_model, direction = "backward", trace = TRUE)
summary(backward_model)

# Cox-PH assumption test
uji.ph <- cox.zph(forward_model)
uji.ph
ggcoxzph(uji.ph)

# Log-log plot
loglog_vars <- c("ulcer", "sex", "age_category", "thickness_strata")

lapply(loglog_vars, function(var) {
  plot(survfit(as.formula(paste("Surv(time, status) ~", var)), data = data), 
       fun = "cloglog", lty = 1.2, mark.time = FALSE, 
       xlab = paste("Survival Time for", gsub("_", " ", var)), 
       ylab = "log(H(t))", xlim = c(1000, 2000))
})

# Partial Likelihood
model_breslow <- coxph(Surv(time, status) ~ ulcer + sex + age_category + thickness_strata, data = data, ties = "breslow")
model_breslow$coefficients

donly <- subset(data, time >= 60 & status == 1)
donly

tht_i <- exp(predict(model_breslow, donly[1,]))
tht_j <- by(donly, seq_len(nrow(donly)), function(row) exp(predict(model_breslow, row)))

tht_i
sum(tht_j)

ptl <- tht_i / sum(tht_j)
ptl
