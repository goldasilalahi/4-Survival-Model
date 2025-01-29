# 1. Import Library ------------------------------------------------------------
library(survival)
library(KMsurv)
library(MASS)
library(survminer)
library(readxl)
library(ggplot2)
library(ggfortify)
library(dplyr)  # untuk manipulasi data
library(car)

# 2. Load Data -----------------------------------------------------------------
data <- read_excel("Library/Containers/com.apple.AppleMediaServicesUI.SpyglassPurchases/Data/file:/Users/Golda/College/sem 4/Model Survival/final project/melanoma.xlsx")

# 3. Remove rows with status 3 -------------------------------------------------
data <- data %>% filter(status != 3)

head(data)
str(data)
dim(data)
summary(data)

# 4. Mendefinisikan variabel faktor --------------------------------------------
data$sex <- as.factor(data$sex)
data$ulcer <- as.factor(data$ulcer)

# 5. Mengubah variabel numerik menjadi variabel kategorik (age dan thickness) --
data$thickness_strata <- cut(data$thickness, 
                             breaks = c(-Inf, 1, 2, 4, Inf), 
                             labels = c("Level 1", "Level 2", "Level 3", "Level 4"),
                             right = FALSE)

data$age_category <- cut(data$age, 
                         breaks = c(4, 11, 25, 45, 65, Inf), 
                         labels = c("kanak-kanak", "remaja", "dewasa", "lansia", "manula"),
                         right = FALSE)

# 6. Plot variabel ---------------------------------------------------------
plots <- list()

# Plot untuk time
plots[['time']] <- ggplot(data, aes(x = time)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black") +
  ggtitle("Distribution of Time")
plots[['time']]

# Plot untuk status
plots[['status']] <- ggplot(data, aes(x = factor(status))) +
  geom_bar(fill = "lightcoral", color = "black") +
  ggtitle("Distribution of Status")
plots[['status']]

# Plot untuk sex
plots[['sex']] <- ggplot(data, aes(x = factor(sex))) +
  geom_bar(fill = "maroon", color = "black") +
  ggtitle("Distribution of Sex")
plots[['sex']]

# Plot untuk age
plots[['age_category']] <- ggplot(data, aes(x = age_category)) +
  geom_bar(fill = "mediumseagreen", color = "black") +
  ggtitle("Distribution of Age")
plots[['age_category']]

# Plot untuk year
plots[['year']] <- ggplot(data, aes(x = year)) +
  geom_histogram(binwidth = 1, fill = "gold", color = "black") +
  ggtitle("Distribution of Year")
plots[['year']]

# Plot untuk thickness
plots[['thickness_strata']] <- ggplot(data, aes(x = thickness_strata)) +
  geom_bar(fill = "turquoise", color = "black") +
  ggtitle("Distribution of Thickness")
plots[['thickness_strata']]

# Plot untuk ulcer
plots[['ulcer']] <- ggplot(data, aes(x = factor(ulcer))) +
  geom_bar(fill = "lightsalmon", color = "black") +
  ggtitle("Distribution of Ulcer")
plots[['ulcer']]

# Tampilkan semua plot
library(gridExtra)
grid.arrange(grobs = plots, ncol = 2)

# 7. Plot Kaplan-Meier -----------------------------------------------------------
# Ulcer
St.ulcer <- surv_fit(Surv(time, status == 1) ~ ulcer, data = data)
ggsurvplot(St.ulcer, data = data, conf.int = TRUE, pval = TRUE,
           risk.table = TRUE,
           title = "", risk.table.height = .3)
# Sex
St.sex <- surv_fit(Surv(time, status == 1) ~ sex, data = data)
ggsurvplot(St.sex, data = data, conf.int = TRUE, pval = TRUE,
           risk.table = TRUE,
           title = "", risk.table.height = .3)
# Age_Category
St.age_category <- surv_fit(Surv(time, status == 1) ~ age_category, data = data)
ggsurvplot(St.age_category, data = data, conf.int = TRUE, pval = TRUE,
           risk.table = TRUE,
           title = "", risk.table.height = .3)
# Thickness_Strata
St.thickness_strata <- surv_fit(Surv(time, status == 1) ~ thickness_strata, data = data)
ggsurvplot(St.thickness_strata, data = data, conf.int = TRUE, pval = TRUE,
           risk.table = TRUE,
           title = "", risk.table.height = .3)

# 8. Log-rank ---------------------------------------------------------------------
surv_object <- Surv(time = data$time, event = data$status == 1)

logrank_sex <- survdiff(surv_object ~ sex, data = data)
logrank_ulcer <- survdiff(surv_object ~ ulcer, data = data)
logrank_thickness_strata <- survdiff(surv_object ~ thickness_strata, data = data)
logrank_age_category <- survdiff(surv_object ~ age_category, data = data)

logrank_sex
logrank_ulcer
logrank_thickness_strata
logrank_age_category

# 9. Analisis Cox-Ph
model_ulcer <- coxph(Surv(time, status) ~ ulcer, data = data)
model_sex <- coxph(Surv(time, status) ~ sex, data = data)
model_age_category <- coxph(Surv(time, status) ~ age_category, data = data)
model_thickness_strata <- coxph(Surv(time, status) ~ thickness_strata, data = data)

model_ulcer
model_sex
model_age_category
model_thickness_strata

ggforest(model_ulcer, data=data)
ggforest(model_sex, data=data)
ggforest(model_age_category, data=data)
ggforest(model_thickness_strata,data=data)

# 10. Mencari kandidat model -------------------------------------------------------
# Membangun model penuh menggunakan Cox Proportional Hazards
full_model <- coxph(surv_object ~ sex + age_category + year + thickness_strata + ulcer, data = data)
summary(full_model)

both_model <- stepAIC(full_model, direction = "both", trace = TRUE)
summary(both_model)

forward_model <- stepAIC(full_model, direction = "forward", trace = TRUE)
summary(forward_model)

backward_model <- stepAIC(full_model, direction = "backward", trace = TRUE)
summary(backward_model)

# 11. Uji asumsi cox-ph ------------------------------------------------------------
# Dengan Residual Schoenfeld 
uji.ph <- cox.zph(forward_model)
uji.ph
ggcoxzph(uji.ph)

# Dengan Plot Log-Log
plot(survfit(Surv(time,status) ~ ulcer, data = data),fun = "cloglog",lty = 1.2,
     mark.time = FALSE, xlabs = "Waktu Survival Ulcer, T", ylab = "log(H(t))", xlim = c(1000, 2000))
plot(survfit(Surv(time,status) ~ sex, data = data),fun = "cloglog",lty = 1.2,
     mark.time = FALSE, xlabs = "Waktu Survival Sex, T", ylab = "log(H(t))", xlim = c(1000, 2000))
plot(survfit(Surv(time,status) ~ age_category, data = data),fun = "cloglog",lty = 1.2,
     mark.time = FALSE, xlabs = "Waktu Survival Age_Category, T", ylab = "log(H(t))", xlim = c(1000, 2000))
plot(survfit(Surv(time,status) ~ thickness_strata, data = data),fun = "cloglog",lty = 1.2,
     mark.time = FALSE, xlabs = "Waktu Survival Thickness_Strata, T", ylab = "log(H(t))", xlim = c(1000, 2000))

# 12. Analisis Partial Likelihood ---------------------------------------------
#Partial Likelihood pada satu titik waktu (non ties)
model_breslow <- coxph(Surv(time, status) ~ ulcer + sex + age_category + 
                       thickness_strata, data = data, ties = "breslow")
model_breslow$coefficients
donly<-subset(data, ...1 >= 60 & status == 1)
donly

tht_i<-exp(predict(model_breslow, donly[1,]))
tht_j<-by(donly, seq_len(nrow(donly)), function(row) 
  exp(predict(model_breslow, row)))
tht_i
sum(tht_j)

ptl<-tht_i/sum(tht_j)
ptl