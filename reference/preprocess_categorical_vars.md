# Preprocess categorical variables (character/factor) to integer codes and dummies

Preprocess categorical variables (character/factor) to integer codes and
dummies

## Usage

``` r
preprocess_categorical_vars(
  data,
  target_vars = NULL,
  dummy_vars = NULL,
  exclude_cols = NULL,
  quiet = FALSE,
  expand_ordered = FALSE
)
```

## Arguments

- data:

  A data.frame or list of data.frames

- quiet:

  Logical; whether to suppress informational messages

## Value

The modified data object with categorical_vars attribute
