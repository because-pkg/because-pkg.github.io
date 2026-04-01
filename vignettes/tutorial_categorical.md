# Tutorial: Handling Categorical & Ordinal Data in `because`

This guide explains how the `because` package processes non-continuous data, how to set up your models, and how to interpret the resulting marginal effects and DAG visualizations.

---

## 1. Data Preparation in R
The behavior of `because` depends heavily on how you define your variables in your R `data.frame`.

| R Data Type | `because` Interpretation | Visualization Behavior |
| :--- | :--- | :--- |
| `factor(x)` | **Unordered Categorical** (Multinomial) | Shows **Bundles** (one arrow per category if `multinomial_probabilities=TRUE`). |
| `factor(x, ordered=T)` | **Ordered Categorical** (Ordinal) | Shows a **Single Arrow** representing a "one-step" increase. |
| `numeric(x)` | **Continuous** (Gaussian) | Shows a **Single Arrow** representing a +1 unit change. |

### Example Setup
```r
# Unordered: Each group is a separate category (D-separation uses dummy vars)
dat$IncomeActivity <- factor(dat$IncomeActivity) 

# Ordered: Level 1 < Level 2 < Level 3 (Uses polynomial contrasts)
dat$Education <- factor(dat$Education, ordered = TRUE)
```

---

## 2. Setting up the `because()` Model
When you fit the model, you must specify the `family` for your **dependent** variables (Responses).

```r
fit <- because(
  equations = list(
    Support ~ Education,
    Education ~ IncomeActivity
  ),
  data = dat,
  family = c(
    Support = "binomial",    # Logistic regression (0/1)
    Education = "ordinal"    # Proportional odds model (1..K)
  )
)
```

### How `because` behaves:
1.  **As a Response**: If a variable is in `family = "ordinal"`, `because` automatically uses a **Latent Variable approach**. It assumes there is a hidden continuous scale behind the categories and estimates "thresholds" (cut-points) between the levels.
2.  **As a Predictor**: 
    *   **If Ordinal**: It uses **Orthogonal Polynomials** (Linear, Quadratic, etc.). The "Marginal Effect" is the average change in the outcome when moving from Category $k$ to $k+1$.
    *   **If Unordered Factor**: It uses **Dummy Coding**. Every category is compared against the first one (Reference).

---

## 3. Average Marginal Effects (AME)
The `marginal_effects(fit)` function converts complex non-linear parameters into an intuitive scale: **the expected change in the probability/value of the outcome.**

> [!IMPORTANT]
> **Deterministic Consistency**: Because marginal effects are calculated post-hoc using posterior samples, `because` uses a **fixed internal seed** and a high-precision default (**500 samples**). This ensures that your DAG plot and your Coefficient plot always show identical numerical results.

### The Calculation Logic:
1.  **Base Prediction**: It calculates the current expected outcome for the whole dataset.
2.  **The Intervention**: It manually "shifts" the predictor (e.g., everyone moves up one education level).
3.  **The Difference**: It recalculates the expected outcome and takes the average difference (`New` - `Base`).

---

## 4. Comparing Visualization Tools

| Feature | `plot_dag()` | `plot_coef()` |
| :--- | :--- | :--- |
| **Purpose** | Shows the "Big Picture" structure and flow. | Shows "Publication-Ready" precision for every path. |
| **Multinomials** | Uses **Bundles** (Arc groups) for categories. | Uses **Individual Rows** (One per category). |
| **Significance** | Arrow thickness + color (Black = Sig). | Point-and-whisker (Black = Sig). |
| **Default Scale** | Marginal Effects (AME). | Supports both Marginal and Raw Beta. |

### Why an effect might look different:
*   **Marginal vs Raw**: A raw coefficient ($\beta$) might have a wide interval, but its **Marginal Effect** (the actual impact on the outcome) can be significant because it accounts for the distribution of other variables in your dataset.
*   **Significance**: In both plots, **Black signifies significance** (the 95% Credible Interval does not cross zero). Grey signifies a non-significant effect.

---

## 5. Global Fit (d-separation)
The `dsep` tests in `because` are computationally rigorous. They do **not** assume your data is linear.
*   The package preserves your `categorical_vars` metadata during these tests, ensuring that even in small sub-models, the number of levels ($K$) and threshold structures are correctly maintained.

> [!TIP]
> **Pro Tip**: To generate perfectly synced outputs for your manuscript, use:
> ```r
> p1 <- plot_dag(fit, type = "marginal")
> p2 <- plot_coef(fit, type = "marginal")
> ```
