# Prepare Multiscale Data for JAGS

Transforms multiscale dataframes into a flat list of vectors and ID
indices suitable for JAGS

## Usage

``` r
prepare_hierarchical_jags_data(hierarchical_info, vars_needed)
```

## Arguments

- hierarchical_info:

  List containing 'data', 'levels', 'hierarchy', 'link_vars'

- vars_needed:

  Character vector of all variable names needed for the model

## Value

List with 'data_list' (for jags) and 'n_vec' (named vector of sample
sizes)
