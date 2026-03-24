# Monitor Variables Hook

Modules implement this to specify which parameters to monitor for a
given response (e.g. z_Y for occupancy).

## Usage

``` r
get_monitor_vars_hook(family, variable_name, ...)
```

## Arguments

- family:

  The S3 family list.

- variable_name:

  Name of the response variable.

- ...:

  Additional arguments.

## Value

Character vector of parameter names to monitor.
