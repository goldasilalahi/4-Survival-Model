# 4-Survival-Model
This project applies **survival analysis** techniques to melanoma patient data. The methods used include:  
- **Kaplan-Meier Estimation**: Non-parametric method to estimate survival probabilities.  
- **Log-Rank Test**: Hypothesis test to compare survival distributions across groups.  
- **Cox Proportional Hazards Model**: Regression model to estimate hazard ratios and identify significant risk factors.  
- **Model Selection (AIC-based Stepwise Selection)**: Determines the best-fitting Cox model using forward, backward, and stepwise approaches.  
- **Proportional Hazards Assumption Testing**: Validates the Cox model using Schoenfeld residuals and log-log plots.  

## Output
- **Kaplan-Meier curves** show lower survival for patients with ulcers and higher tumor thickness.  
- **Log-rank test** confirms significant differences in survival for ulcer presence (**p < 0.001**) and tumor thickness levels.  
- **Cox model results** indicate that **tumor thickness and ulcer presence significantly impact survival**.  
- **Final model selection** (stepwise AIC) retains **tumor thickness and ulcer presence as key predictors**.
