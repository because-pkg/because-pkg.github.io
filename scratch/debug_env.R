# Debug Generic Visibility
library(because)

cat("--- Debugging Generic Visibility ---\n")
gen_name <- "jags_structure_definition"

# 1. Check Search Path
cat("Exists in Search Path? ", exists(gen_name), "\n")
if (exists(gen_name)) {
  obj <- get(gen_name)
  cat("Environment of Object: ", capture.output(print(environment(obj))), "\n")
}

# 2. Check Namespace
ns <- asNamespace("because")
cat("Exists in 'because' Namespace? ", exists(gen_name, envir = ns), "\n")
if (exists(gen_name, envir = ns)) {
  cat("Object in Namespace found.\n")
}

# 3. Check Exports
cat("Is among Namespace Exports? ", gen_name %in% getNamespaceExports(ns), "\n")

# 4. Try manual registration test
test_method <- function(...) cat("Test success\n")
cat("Attempting manual registerS3method test...\n")
tryCatch({
  registerS3method(gen_name, "debug_class", test_method, envir = ns)
  cat("Registration test into 'because' namespace: SUCCESS\n")
}, error = function(e) {
  cat("Registration test into 'because' namespace: FAILED - ", e$message, "\n")
})

cat("--- Debug finished ---\n")
