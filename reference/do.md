# Causal Interventions for Structural Equation Models

Evaluates counterfactuals using Pearl's do-operator. Supports atomic
interventions (setting a variable to a fixed value), shift interventions
(adding a constant to the natural value), and stochastic interventions
(drawing from a distribution).

## Usage

``` r
do(object, ...)

because_do(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Interventions. Can be numeric constants, formulas, or expressions.

## Value

Intervened data or counterfactual samples.
