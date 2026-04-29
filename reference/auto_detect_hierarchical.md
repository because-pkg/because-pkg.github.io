# Auto-Detect Multiscale Data Structure

Infers scales, hierarchy, and link variables from a list of data.frames

## Usage

``` r
auto_detect_hierarchical(data, eq_vars, quiet = FALSE)
```

## Arguments

- data:

  List of data.frames at different scales

- eq_vars:

  Character vector of variable names used in equations

- quiet:

  Logical, if TRUE suppress messages

## Value

List with 'levels', 'hierarchy', 'link_vars'
