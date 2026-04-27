if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "Estimate",
        "LowerCI",
        "Test",
        "UpperCI",
        "edge_id",
        "edge_label",
        "edge_type",
        "label_display",
        "model_label",
        "occ_species",
        "significant",
        "weight_abs",
        "xend",
        "xmax",
        "xmin",
        "y",
        "yend",
        "ymax",
        "ymin",
        "Estimate",
        "LowerCI",
        "Test",
        "UpperCI", # repeated for clarity if needed
        "name",
        "name",
        "type",
        "Label",
        "Lower",
        "Path",
        "Significant",
        "Upper",
        "curvature",
        "dx",
        "dy",
        "to",
        "returnType"
    ))
}

# Internal environment for storing dynamic functions (e.g. for NIMBLE)
# to avoid polluting the global environment and satisfy CRAN requirements.
.because_env <- new.env(parent = emptyenv())
