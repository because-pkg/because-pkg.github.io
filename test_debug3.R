library(because)
library(because.phybase)
data(rhino.dat)
data(rhino.tree)
sem8_eq <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

# Let's call the internal function to see what it returns
cat("Testing is_valid_random_level:\n")
# because:::is_valid_random_level is not exported.
h_info <- because:::process_hierarchical_data(rhino.dat, id_col="SP")
cat("LS:", because:::is_valid_random_level("LS", "SP", h_info), "\n")
cat("BM:", because:::is_valid_random_level("BM", "SP", h_info), "\n")
