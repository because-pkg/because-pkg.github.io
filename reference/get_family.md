# Create a family object for a given distribution

Returns a family object that can be dispatched on via S3.

## Usage

``` r
get_family(name)
```

## Arguments

- name:

  Character string; family name (e.g., "gaussian", "binomial")

## Value

A family object with appropriate S3 class

## Examples

``` r
get_family("gaussian")
#> $family
#> [1] "gaussian"
#> 
#> attr(,"class")
#> [1] "because_family_gaussian" "because_family"         
get_family("binomial")
#> $family
#> [1] "binomial"
#> 
#> attr(,"class")
#> [1] "because_family_binomial" "because_family"         
```
