# vi: fdm=marker
# Script needed to run testthat automatically from ‘R CMD check’. See testthat::test_dir documentation.
library(testthat)
library(mineMS2)

# Run tests {{{1
################################################################

Sys.setenv(TESTTHAT_REPORTER="summary")
test_check("mineMS2")
