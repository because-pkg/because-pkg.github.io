library(because)

cat("=== Testing Polynomial Term Parsing ===\n\n")

# Test 1: Simple polynomial
cat("Test 1: Simple polynomial\n")
f1 <- y ~ x + I(x^2)
poly1 <- because:::extract_polynomial_terms(f1)
print(poly1)
cat("\n")

# Test 2: Multiple powers
cat("Test 2: Multiple powers\n")
f2 <- y ~ x + I(x^2) + I(x^3)
poly2 <- because:::extract_polynomial_terms(f2)
print(poly2)
cat("\n")

# Test 3: Multiple variables
cat("Test 3: Multiple variables\n")
f3 <- y ~ x + I(x^2) + z + I(z^2)
poly3 <- because:::extract_polynomial_terms(f3)
print(poly3)
cat("\n")

# Test 4: Formula expansion
cat("Test 4: Formula expansion\n")
f4 <- y ~ x + I(x^2) + I(x^3)
poly4 <- because:::extract_polynomial_terms(f4)
expanded <- because:::expand_polynomial_formula(f4, poly4)
cat("Original: ", deparse(f4), "\n")
cat("Expanded: ", deparse(expanded), "\n")
cat("\n")

# Test 5: JAGS generation
cat("Test 5: JAGS code generation\n")
jags_code <- because:::generate_polynomial_jags(poly4)
cat(paste(jags_code, collapse = "\n"), "\n")
cat("\n")

# Test 6: No polynomials
cat("Test 6: No polynomials (should return NULL)\n")
f5 <- y ~ x + z
poly5 <- because:::extract_polynomial_terms(f5)
print(poly5)
cat("\n")

cat("=== All Parsing Tests Complete! ===\n")
