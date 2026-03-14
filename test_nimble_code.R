library(nimble)
model_string <- "
  for(i in 1:10) {
    y[i] ~ dnorm(mu, 1)
  }
  mu ~ dnorm(0, 1)
"

# Method: nimbleCode requires an expression, parse creates an expression list
# So we wrap the string in '{' '}'
safe_string <- paste0("{", model_string, "}")
code <- parse(text = safe_string)[[1]]

print("Code generated:")
print(code)

tryCatch(
  {
    m1 <- nimbleModel(code = code)
    print("Method 1 succeeded")
  },
  error = function(e) print(paste("Method 1 failed:", e$message))
)
