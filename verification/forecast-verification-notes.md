# Forecast Verification: Methods and Key Concepts

> **Source:** [CAWCR — WWRP/WGNE Joint Working Group on Forecast Verification Research](https://www.cawcr.gov.au/projects/verification/verif_web_page.html)

---

## Table of Contents

1. [What is Forecast Verification?](#1-what-is-forecast-verification)
2. [Why Verify Forecasts?](#2-why-verify-forecasts)
3. [Types of Forecasts](#3-types-of-forecasts)
4. [What Makes a Forecast "Good"?](#4-what-makes-a-forecast-good)
5. [Forecast Quality vs. Value](#5-forecast-quality-vs-value)
6. [Standard Verification Methods](#6-standard-verification-methods)
   - [Dichotomous (Yes/No) Forecasts](#61-methods-for-dichotomous-yesno-forecasts)
   - [Multi-Category Forecasts](#62-methods-for-multi-category-forecasts)
   - [Continuous Forecasts](#63-methods-for-continuous-forecasts)
   - [Probabilistic Forecasts](#64-methods-for-probabilistic-forecasts)
7. [Advanced / Diagnostic Methods](#7-advanced--diagnostic-verification-methods)
   - [Spatial Verification](#71-spatial-forecast-verification)
   - [Ensemble / Probabilistic Diagnostics](#72-ensembleprobabilistic-diagnostics)
   - [Rare Event Verification](#73-rare-event-verification)
   - [Other Diagnostic Tools](#74-other-diagnostic-tools)
8. [Practical Considerations](#8-practical-considerations)
9. [Key Takeaways](#9-key-takeaways)

---

## 1. What is Forecast Verification?

Forecast verification is the process of **comparing forecasts against corresponding observations** (the "truth") to assess forecast quality. It can be qualitative (subjective, visual) or quantitative (using statistical metrics).

---

## 2. Why Verify Forecasts?

Three primary reasons:

| Purpose | Description |
|---|---|
| **Administrative/Monitoring** | Track forecast quality over time; detect improvements or degradation |
| **Scientific/Improvement** | Identify strengths and weaknesses to improve forecast systems |
| **Economic/Comparative** | Compare competing forecast systems or methods to pick the best |

---

## 3. Types of Forecasts

Forecasts can be classified along multiple dimensions:

### By Nature
- **Deterministic** — single-valued prediction (e.g., "tomorrow's high will be 30 C")
- **Probabilistic** — probability distribution (e.g., "70% chance of rain")
- **Qualitative** — categorical/descriptive (e.g., "partly cloudy")

### By Space-Time Domain
- **Time series** — sequential forecasts at a fixed point
- **Spatial distribution** — forecast fields over a region (maps)
- **Pooled** — aggregated over time and/or space

### By Specificity
| Type | Description | Example |
|---|---|---|
| **Dichotomous** | Yes/No forecasts | "Will it rain?" |
| **Multi-category** | Several discrete categories | "Light / Moderate / Heavy rain" |
| **Continuous** | Numerical values | Temperature = 25.3 C |
| **Object-oriented** | Spatial features/entities | Storm cells, rainfall areas |

![Time series example](images/timeseries.gif)

![Spatial forecast maps](images/DWDmaps.gif)

---

## 4. What Makes a Forecast "Good"?

Murphy (1993) defined three types of forecast "goodness":

1. **Consistency** — forecasts match the forecaster's best judgment (internal)
2. **Quality** — forecasts match what actually happened (external)
3. **Value** — forecasts help users make better decisions (economic)

### Nine Attributes of Forecast Quality

| Attribute | What It Measures |
|---|---|
| **Bias** | Systematic tendency to over-forecast or under-forecast |
| **Association** | Strength of linear relationship between forecasts and observations |
| **Accuracy** | Average correspondence between forecasts and observations |
| **Skill** | Accuracy relative to a reference (e.g., climatology or persistence) |
| **Reliability** | Whether stated probabilities match observed frequencies |
| **Resolution** | Ability to distinguish different observed outcomes |
| **Sharpness** | Tendency to forecast extreme (confident) probabilities |
| **Discrimination** | Ability to distinguish between event occurrences and non-occurrences |
| **Uncertainty** | Variability of the observations themselves |

---

## 5. Forecast Quality vs. Value

- **Quality** = how well forecasts correspond to observations (measured by statistical metrics)
- **Value** = whether forecasts help decision-makers achieve better outcomes (measured in economic terms)

> A high-quality forecast may have **no value** if the user cannot act on it.
> Conversely, a lower-quality forecast can still have **high value** if it informs a critical decision.

![Relative value — deterministic](images/value_determ.gif)

![Relative value — ensemble](images/value_EPS.gif)

---

## 6. Standard Verification Methods

### 6.0 Visual (Eyeball) Verification

Examining forecasts and observations side by side using human judgment. Intuitive and flexible, but **subjective and prone to bias**.

---

### 6.1 Methods for Dichotomous (Yes/No) Forecasts

#### The Contingency Table

The foundation of dichotomous verification:

|  | **Observed: Yes** | **Observed: No** |
|---|---|---|
| **Forecast: Yes** | **a** (Hits) | **b** (False Alarms) |
| **Forecast: No** | **c** (Misses) | **d** (Correct Negatives) |

Total: `n = a + b + c + d`

#### Key Metrics

| Metric | Formula | Range | Perfect | Notes |
|---|---|---|---|---|
| **Accuracy (Fraction Correct)** | `(a + d) / n` | 0 to 1 | 1 | Misleading for rare events |
| **Bias Score** | `(a + b) / (a + c)` | 0 to inf | 1 | >1 = overforecast, <1 = underforecast |
| **Probability of Detection (POD)** | `a / (a + c)` | 0 to 1 | 1 | Sensitivity; ignores false alarms |
| **False Alarm Ratio (FAR)** | `b / (a + b)` | 0 to 1 | 0 | Fraction of "yes" forecasts that were wrong |
| **Prob. of False Detection (POFD)** | `b / (b + d)` | 0 to 1 | 0 | False alarm rate among non-events |
| **Threat Score (CSI)** | `a / (a + b + c)` | 0 to 1 | 1 | Ignores correct negatives; good for rare events |
| **Equitable Threat Score (ETS)** | `(a - a_ref) / (a + b + c - a_ref)` | -1/3 to 1 | 1 | Adjusted for random chance hits |
| **Hanssen-Kuipers (HK)** | `POD - POFD` | -1 to 1 | 1 | a.k.a. True Skill Statistic (TSS), Peirce Skill Score |
| **Heidke Skill Score (HSS)** | `2(ad - bc) / [(a+c)(c+d) + (a+b)(b+d)]` | -1 to 1 | 1 | Skill relative to random chance |
| **Odds Ratio** | `(a * d) / (b * c)` | 0 to inf | inf | Ratio of odds of a hit to odds of a false alarm |
| **Odds Ratio Skill Score (Yule's Q)** | `(ad - bc) / (ad + bc)` | -1 to 1 | 1 | Base-rate independent |

Where `a_ref = (a + b)(a + c) / n` (expected random hits)

#### The Finley Tornado Example (1884)

A classic case study illustrating the **pitfall of using accuracy alone**:

| | Observed: Tornado | Observed: No Tornado |
|---|---|---|
| Forecast: Tornado | 28 (hits) | 72 (false alarms) |
| Forecast: No Tornado | 23 (misses) | 2680 (correct negatives) |

- **Accuracy = 96.6%** — sounds great!
- But a "never forecast tornado" strategy yields **96.4%** accuracy
- The forecast adds almost no skill despite high accuracy
- **Lesson:** For rare events, accuracy is misleading. Use POD, FAR, CSI, ETS, or skill scores instead.

---

### 6.2 Methods for Multi-Category Forecasts

Extension of the contingency table to multiple categories (e.g., light/moderate/heavy rain).

![Category histogram](images/histogram.gif)

#### Key Approaches

- **Distributions approach** — examine the full contingency table; perfect forecasts have counts only on the diagonal
- **Histogram comparison** — compare forecast vs. observed frequency distributions
- **Multi-category Accuracy** — fraction of forecasts on the diagonal
- **Heidke Skill Score (multi-cat)** — accuracy relative to random chance
- **Hanssen-Kuipers (multi-cat)** — uses unbiased random baseline

---

### 6.3 Methods for Continuous Forecasts

For forecasts of numerical values (temperature, wind speed, pressure, etc.)

#### Exploratory Plots

![Scatter plot](images/scatterplot.gif)

![Error scatter plot](images/error_scatterplot.gif)

![Box plot](images/boxplot.gif)

#### Summary Metrics

| Metric | Formula | Range | Perfect | Notes |
|---|---|---|---|---|
| **Mean Error (ME / Bias)** | `mean(F - O)` | -inf to +inf | 0 | Additive bias; errors cancel out |
| **Multiplicative Bias** | `mean(F) / mean(O)` | 0 to inf | 1 | Best for bounded-at-zero quantities |
| **Mean Absolute Error (MAE)** | `mean(\|F - O\|)` | 0 to inf | 0 | Average error magnitude |
| **Root Mean Square Error (RMSE)** | `sqrt(mean((F-O)^2))` | 0 to inf | 0 | Penalizes large errors more |
| **Correlation Coefficient (r)** | Pearson's r | -1 to 1 | 1 | Linear association; ignores bias |
| **Anomaly Correlation** | r computed on anomalies from climatology | -1 to 1 | 1 | Removes climatological mean |
| **LEPS** | Error in cumulative probability space | 0 to 1 | 0 | Accounts for climatological distribution |
| **S1 Score** | Gradient error measure | 0 to 200 | 0 | Resolution-dependent |
| **Skill Score** | `1 - (MSE_fcst / MSE_ref)` | -inf to 1 | 1 | Improvement over reference forecast |

#### MSE Decomposition

MSE can be decomposed into interpretable components:

```
MSE = Bias^2 + Variance_error
    = (mean(F) - mean(O))^2 + (std(F) - std(O))^2 + 2*std(F)*std(O)*(1 - r)
```

This separates error into **unconditional bias**, **conditional bias** (amplitude error), and **phase/pattern error**.

---

### 6.4 Methods for Probabilistic Forecasts

Probabilistic forecasts express uncertainty as probabilities (e.g., 70% chance of rain).

#### Three Desirable Properties

1. **Reliability** — forecast probabilities match observed frequencies (e.g., of all 70% rain forecasts, ~70% should verify)
2. **Sharpness** — probabilities are near 0 or 1 (confident), not clustered near climatology
3. **Resolution** — the system distinguishes occasions when the event occurs from when it does not

#### Reliability Diagram

![Reliability diagram](images/ReliabilityDiagram.gif)

- X-axis: forecast probability bins
- Y-axis: observed relative frequency
- **Perfect reliability** = points on the diagonal
- Above diagonal = under-forecasting
- Below diagonal = over-forecasting

#### Brier Score (BS)

```
BS = (1/N) * SUM(f_i - o_i)^2
```

Where `f_i` = forecast probability, `o_i` = observation (0 or 1)

- Range: 0 to 1
- Perfect: 0
- Can be decomposed: `BS = Reliability - Resolution + Uncertainty`

#### Brier Skill Score (BSS)

```
BSS = 1 - (BS / BS_climatology)
```

- Range: -inf to 1
- Perfect: 1; Zero = no skill over climatology

#### ROC (Relative Operating Characteristic)

![ROC curve](images/ROC.gif)

- Plots **Hit Rate (POD)** vs. **False Alarm Rate (POFD)** at varying probability thresholds
- **Area under curve (AUC):** 1.0 = perfect discrimination, 0.5 = no skill (diagonal)
- Measures **discrimination** (not reliability)

#### Ranked Probability Score (RPS)

For multi-category probabilistic forecasts. Penalizes predictions further from the actual outcome category.

```
RPS = (1/J) * SUM_j (CDF_forecast_j - CDF_observed_j)^2
```

#### Relative Value Score

Measures economic benefit for a given cost-to-loss ratio. Plotted as a function of the decision threshold.

---

## 7. Advanced / Diagnostic Verification Methods

### 7.1 Spatial Forecast Verification

Traditional point-by-point verification suffers from the **"double penalty" problem** — a displaced forecast is penalized twice (miss at correct location + false alarm at forecast location). Spatial methods address this.

#### Scale Decomposition

Uses **wavelet** or **discrete cosine transforms** to separate forecast errors by spatial scale. Reveals whether errors are at large scales (synoptic) or small scales (mesoscale).

![Intensity-scale MSE skill](images/IS_MSEskill.gif)

#### Fuzzy / Neighborhood Methods

Relax the exact-match requirement by considering neighborhoods around grid points.

- **Fractions Skill Score (FSS):** Compares fractional coverage within expanding windows
  - FSS = 0: no skill
  - FSS = 1: perfect
  - Identifies the **minimum useful scale** of a forecast

![Fractions Skill Score plot](images/FSSplot.gif)

#### Object-Oriented Methods

Identify discrete features (e.g., storm cells, rainfall blobs) and verify their properties.

**CRA (Contiguous Rain Area):**
- Decomposes total error into **location error**, **volume error**, and **pattern error**

![CRA entities](images/CRA_entities0.gif)

**MODE (Method for Object-based Diagnostic Evaluation):**
- Identifies objects in both forecast and observation fields
- Matches objects based on proximity, overlap, intensity
- Evaluates attributes: location, size, shape, intensity, orientation

![MODE objects](images/MODE.gif)

**SAL (Structure-Amplitude-Location):**
- Three independent components each ranging from -2 to +2
- S = Structure, A = Amplitude, L = Location
- Perfect: S = A = L = 0

---

### 7.2 Ensemble/Probabilistic Diagnostics

#### Rank Histogram (Talagrand Diagram)

![Rank histogram](images/RankHistogram.gif)

Evaluates whether ensemble spread is adequate by plotting where the observation falls within the ranked ensemble members.

| Shape | Interpretation |
|---|---|
| **Flat** | Ensemble spread is appropriate |
| **U-shaped** | Ensemble spread is too narrow (underdispersive) |
| **Dome-shaped** | Ensemble spread is too wide (overdispersive) |
| **Asymmetric** | Ensemble has a systematic bias |

#### Correspondence Ratio

![Correspondence ratio](images/CorrespondenceRatio.gif)

![Venn diagram](images/Venn.gif)

Ratio of intersection to union area for ensemble member agreement.

#### Ignorance Score (Logarithmic Scoring Rule)

```
Ignorance = -log2(p_observed)
```

Where `p_observed` = forecast probability assigned to the outcome that actually occurred. Heavily penalizes confident wrong predictions.

---

### 7.3 Rare Event Verification

Rare events pose special challenges — standard metrics become unreliable with small sample sizes.

![Deterministic limit](images/DetLim.gif)

- **Deterministic Limit:** Identifies the forecast lead time where skill drops to 0.5
- **Extreme Dependency Score (EDS):** Designed to be stable for rare events; measures association independent of bias
- **Probability Model Approach:** Fits parametric distributions to reduce sampling noise

---

### 7.4 Other Diagnostic Tools

#### Taylor Diagram

Summarizes three statistics on a single polar plot:
- **Correlation coefficient** (angular position)
- **Standard deviation ratio** (radial distance)
- **Centered RMSD** (distance from reference point)

Useful for comparing multiple models at a glance.

#### Root Mean Squared Factor (RMSF)

![RMSF equation](images/RMSF.gif)

A multiplicative error metric using logarithmic transformation. Appropriate when errors scale with magnitude.

#### Quantile-Based Categorical Statistics

![Quantile-based categorical scores](images/fall_bothint_querrclim.gif)

Uses statistical quantiles to define categories, ensuring marginal totals are fixed. Removes the influence of climatological frequency on verification scores.

---

## 8. Practical Considerations

### Sample Size
- Larger samples give more trustworthy verification results
- Confidence intervals are essential, especially for rare events
- Rule of thumb: need many more samples for probabilistic verification than deterministic

### Pooling vs. Stratifying
| Approach | Pros | Cons |
|---|---|---|
| **Pooling** | Statistical reliability; larger samples | Masks performance variations across regimes |
| **Stratifying** | Reveals regime-dependent skill | Requires large datasets; risk of small samples |

Best practice: do both. Pool for overall assessment, stratify for diagnostic insight.

### What is "Truth"?
- Observations are the standard verification baseline, but they have measurement errors
- Point observations vs. gridded analyses create representation mismatches
- Always consider the uncertainty in verification data

### The Double-Penalty Problem
Higher-resolution models often verify **worse** than lower-resolution models on traditional point metrics because:
1. A displaced feature is penalized at **both** locations (miss + false alarm)
2. Lower-resolution (smoother) models get penalized less because they forecast broad averages

**Solution:** Use spatial/neighborhood/object-oriented verification methods.

### Hedging
Some verification metrics can be "gamed" by strategic forecasting (hedging). Desirable properties of metrics:
- **Equitable:** random and constant forecasts receive the same (zero) score
- **Proper:** optimal score is achieved by the forecaster's honest best estimate

---

## 9. Key Takeaways

1. **No single metric tells the whole story** — always use multiple complementary metrics
2. **Match the metric to the forecast type** — dichotomous, multi-category, continuous, or probabilistic
3. **Accuracy is misleading for rare events** — use skill scores, POD, FAR, ETS instead
4. **Distinguish quality from value** — a "good" forecast must be useful, not just accurate
5. **Spatial verification requires specialized methods** — traditional point metrics penalize high-resolution forecasts unfairly
6. **Ensemble spread matters** — use rank histograms and reliability diagrams to diagnose ensemble calibration
7. **Always report confidence intervals** — especially for small samples or rare events
8. **Stratify AND pool** — overall scores hide regime-dependent performance
9. **Consider the "truth"** — observations have errors too; verification is not absolute

---

## References

- Murphy, A.H. (1993). What is a good forecast? An essay on the nature of goodness in weather forecasting. *Weather and Forecasting*, 8, 281-293.
- Finley, J.P. (1884). Tornado predictions. *American Meteorological Journal*, 1, 85-88.
- Wilks, D.S. (2011). *Statistical Methods in the Atmospheric Sciences*. 3rd ed. Academic Press.
- Jolliffe, I.T. and Stephenson, D.B. (2012). *Forecast Verification: A Practitioner's Guide in Atmospheric Science*. 2nd ed. Wiley.
- WWRP/WGNE Joint Working Group on Forecast Verification Research: https://www.cawcr.gov.au/projects/verification/
