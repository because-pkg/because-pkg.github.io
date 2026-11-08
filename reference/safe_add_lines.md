# Add lines to BUGS model while preventing node redefinition

Add lines to BUGS model while preventing node redefinition

## Usage

``` r
safe_add_lines(current_lines, new_lines, declared_nodes = NULL)
```

## Arguments

- current_lines:

  Character vector of current model lines

- new_lines:

  Character vector of new lines to add

- declared_nodes:

  Optional character vector or environment tracking declared node names

## Value

Updated character vector of model lines
