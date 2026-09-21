# Update a Fitted because Model

Re-evaluates a fitted `because` model object with modified arguments
(e.g. new equations, updated data, or changed MCMC settings).

## Usage

``` r
# S3 method for class 'because'
update(object, equations, data, ...)
```

## Arguments

- object:

  A fitted `because` model object.

- equations:

  Optional updated formula list.

- data:

  Optional updated data frame.

- ...:

  Additional arguments to update or override in the original model call.

## Value

A new fitted `because` model object.
