# Get Family Object for S3 Dispatch

Converts a family name string into a family class object for S3
dispatch. This provides a general mechanism for family extensions to be
used.

## Usage

``` r
get_family_object(family_name)
```

## Arguments

- family_name:

  The name of the family (e.g., "occupancy", "gaussian").

## Value

An object of class `because_family_<name>` and `because_family`.
