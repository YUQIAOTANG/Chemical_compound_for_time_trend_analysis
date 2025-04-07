
# Trend Analysis for Chemical Degradation

This project performs kinetic modeling of chemical degradation trends under different environmental conditions, using **natural** and **artificial** trend datasets.

## Input Format

The input file should be an Excel file (`.xlsx`) containing the following columns:

| Name           | Condition        | Concentration | std          | Group         |
|----------------|------------------|----------------|---------------|---------------|
| Name_of_compounds  | 0yr_D_Dry_0d     | 12345     | 67890   | 0yr_D_Dry     |

- `Name`: Name of the chemical compound.
- `Condition`: Encodes time information (e.g., `0yr_D_Dry_0d`, where time is extracted in weeks/days).
- `Concentration`: Observed peak area.
- `std`: Standard deviation of the measurement.
- `Group`: Group label (used for classification and visualization).

---

## Model Overview

Three models are fitted to each compound-group combination, based on the concentration-time profile. These three models are refered from ["Time-concentration profiles of tire particle additives and transformation products under natural and artificial aging"](https://www.sciencedirect.com/science/article/pii/S0048969722072503?via%3Dihub) 

### 1. **Monotonic Decreasing Model**

Used when concentration strictly decreases with time.

**Formula:**

```
C(t) = C₀ * exp(-t / τ)
```

- `C₀`: Initial concentration (fixed from time = 0)
- `τ`: Time constant

---

### 2. **Monotonic Increasing Model**

Used when concentration strictly increases with time.

**Formula:**

```
C(t) = (C₀ - C_f) * exp(-t / τ) + C_f
```

- `C₀`: Initial concentration (fixed from time = 0)
- `C_f`: Final concentration (fixed from max time)
- `τ`: Time constant

---

### 3. **Non-Monotonic Model**

Used when concentration increases and then decreases (or vice versa).

**Formula:**

```
C(t) = C₀ + Cₓ * [-exp(-t / τ₁) + exp(-t / τ₂)]
```

- `C₀`: Initial concentration (fixed)
- `Cₓ`: Amplitude
- `τ₁`, `τ₂`: Time constants

---

## Model Selection Logic

For each group of a compound:

1. **Data Pre-check:**
   - At least 3 time points
   - Must include 0w or 0d measurement

2. **Model Selection:**
   - If all concentration differences ≤ 0 → use **monotonic decreasing**
   - If all ≥ 0 → use **monotonic increasing**
   - Otherwise → try **non-monotonic**

---

## Fallback Conditions (for non-monotonic failure)

If non-monotonic fitting fails (e.g., poor correlation):

- If correlation ≥ 0.7 → fallback to **monotonic increasing**
- If correlation ≤ -0.7 → fallback to **monotonic decreasing**
- Otherwise → log warning, **no model is used**

---

## Output

- Plot `.png` for each compound with error bars, fitted curves, and legends.
- Excel file `kinetic_results.xlsx` with fitted parameters:
  - `Compound`, `Group`, `Model`
  - `C₀`, `Cₓ`, `C_f`, `τ₁`, `τ₂`
  - Standard errors `σ₁`, `σ₂`

---

## Scripts

- `naturaltrend_3.17.R`: Analyzes natural degradation patterns (e.g., "0yr Dry", "0yr Wet").
- `artificialtrend3.24.R`: Analyzes artificial trend under Light/Drak treatments (e.g., "TNU_D_Dry", "TNU_L_Wet").

---

## Dependencies

- `readxl`, `ggplot2`, `dplyr`, `stringr`, `minpack.lm`, `openxlsx`

Install via:

```r
install.packages(c("readxl", "ggplot2", "dplyr", "stringr", "minpack.lm", "openxlsx"))
```

---

## Author

Yuqiao Tang, EAI Northeastern University
