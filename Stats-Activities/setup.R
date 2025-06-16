### Workshop Setup Script ###
### Statistics for Computational Biology Projects Workshop ###
### Run this script before the workshop to install all required packages ###

# Function to install packages if they're not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages) > 0) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, dependencies = TRUE)
  } else {
    cat("All CRAN packages are already installed.\n")
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages) > 0) {
    cat("Installing missing Bioconductor packages:", paste(new_packages, collapse = ", "), "\n")
    BiocManager::install(new_packages, version = "3.18")
  } else {
    cat("All Bioconductor packages are already installed.\n")
  }
}

# CRAN packages needed for the workshop
cran_packages <- c(
  "pwr",           # Power analysis
  "MASS",          # GLM functions
  "ggplot2",       # Data visualization
  "stringr",       # String manipulation
  "lattice",       # Graphics
  "pscl",          # Zero-inflated models
  "lmtest",        # Model testing
  "outliers",      # Outlier detection
  "EnvStats",      # Environmental statistics
  "dbscan"         # Clustering
)

# Bioconductor packages needed for the workshop
bioc_packages <- c(
  "BSgenome",
  "BSgenome.Mmusculus.UCSC.mm9",
  "TxDb.Mmusculus.UCSC.mm9.knownGene",
  "GenVisR"
)

# Install packages
cat("=== Installing CRAN packages ===\n")
install_if_missing(cran_packages)

cat("\n=== Installing Bioconductor packages ===\n")
install_bioc_if_missing(bioc_packages)

# Test that all packages can be loaded
cat("\n=== Testing package loading ===\n")
all_packages <- c(cran_packages, bioc_packages)
failed_packages <- c()

for(pkg in all_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    failed_packages <- c(failed_packages, pkg)
    cat("✗ Failed to load:", pkg, "\n")
  } else {
    cat("✓ Successfully loaded:", pkg, "\n")
  }
}

# Summary
if(length(failed_packages) == 0) {
  cat("\n🎉 All packages installed and loaded successfully!\n")
  cat("You're ready for the workshop.\n")
} else {
  cat("\n⚠️  Some packages failed to install/load:\n")
  cat(paste(failed_packages, collapse = ", "), "\n")
  cat("Please try installing these manually or contact the instructor.\n")
}

# Check R version
cat("\n=== System Information ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")

# Create a session info file for troubleshooting
writeLines(capture.output(sessionInfo()), "session_info.txt")
cat("Session info saved to 'session_info.txt'\n")