# D-Sep Tree Hook

Modules implement this to decide if a tree should be passed to a
specific d-separation test.

## Usage

``` r
dsep_tree_hook(tree, test_eq, hierarchical_info, levels, ...)
```

## Arguments

- tree:

  The structure object (e.g. phylo).

- test_eq:

  The current d-separation test equation.

- hierarchical_info:

  Hierarchical metadata.

- levels:

  Level mapping.

- ...:

  Additional arguments.

## Value

The tree object or NULL.
