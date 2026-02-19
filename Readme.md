# Survival Analysis of Haemodialysis vs Conservative Treatment in ESRD Patients

📊 This project demonstrates advanced survival modeling techniques applied to real-world clinical data.

## Objective

To evaluate the effectiveness of haemodialysis compared to conservative treatment in elderly patients (≥65 years) with end-stage renal disease (ESRD), using advanced survival analysis techniques.

The primary endpoint was time-to-death (all-cause mortality).

---

## Dataset

- N = 8,796 patients
- 174 Conservative treatment
- 8,622 Haemodialysis
- Follow-up up to 16.2 years
- Multiple baseline covariates (age, sex, etc)

---

## Statistical Workflow

### 1. Descriptive Analysis
- Baseline comparison between treatment groups
- Crude cumulative risks
- Mortality rates per person-year
- Incidence Rate Ratio (IRR)

---

### 2. Non-Parametric Survival Analysis
- Kaplan–Meier survival curves
- Log-rank test
- Wilcoxon (Peto–Peto) test
- Conditional survival probabilities

---

### 3. Cox Proportional Hazards Models

#### Univariate Models
- Treatment effect
- Age effect

#### Multivariable Model
Adjusted for:
- Age
- Sex
- COPD
- Diabetes
- Hypertension
- Heart disease
- Liver disease
- Neoplasia
- Vascular disease

Estimated Hazard Ratios (HR) with 95% confidence intervals.

---

### 4. Assessment of PH Assumption

- Log(-log(S(t))) plots
- Schoenfeld residuals
- Grambsch & Therneau test

Result: Evidence of violation of proportional hazards for treatment.

---

### 5. Time-Varying Effects

- Extended Cox model with time-dependent interaction
- Estimated HR(t) and time-varying effectiveness

---

### 6. Piecewise Cox Model

Cut-points:
- 0–6 months
- 6–12 months
- >12 months

Estimated period-specific hazard ratios and effectiveness.

---

### 7. Model Diagnostics

- Cox–Snell residual analysis
- Global goodness-of-fit assessment

---

### 8. Parametric Survival Modeling

- Weibull distribution fit (graphical + MLE)
- Shape and scale parameter estimation
- Multivariable Weibull (AFT form)
- Transformation to PH interpretation

---

## Key Findings

- Crude risk comparison was misleading due to unequal follow-up.
- Person-time analysis showed lower mortality rate under haemodialysis.
- Multivariable Cox model: HR ≈ 0.22 (strong protective effect).
- Proportional hazards assumption violated for treatment.
- Time-varying and piecewise models revealed strong early benefit that attenuates over time.
- Weibull model provided reasonable parametric fit.

---

## Tools Used

- R
- survival package
- survreg

---

## Methods Demonstrated

Kaplan–Meier estimation  
Cox PH modeling  
Time-dependent covariates  
Piecewise PH modeling  
Schoenfeld diagnostics  
Cox–Snell residuals  
Weibull parametric survival modeling  
Hazard ratio interpretation  
Confounding adjustment  

---

## Reproducibility

Code is structured in sequential workflow:
1. Data preparation
2. Descriptive statistics
3. KM analysis
4. Cox modeling
5. Diagnostics
6. Parametric modeling
