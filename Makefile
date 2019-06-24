# vi: fdm=marker

# Global variables {{{1
################################################################

# Zip
PKG_VERSION=$(shell grep '^Version:' DESCRIPTION | sed 's/^Version: //')
ZIPPED_PKG=rhelpers_$(PKG_VERSION).tar.gz

# Set testthat reporter
ifndef TESTTHAT_REPORTER
ifdef VIM
TESTTHAT_REPORTER=summary
else
TESTTHAT_REPORTER=progress
endif
endif

# Default target {{{1
################################################################

all:

# Check and test {{{1
################################################################

check: $(ZIPPED_PKG)
	time R CMD check --no-build-vignettes "$<"

test:
	R -q -e "devtools::test('$(CURDIR)', reporter = c('$(TESTTHAT_REPORTER)', 'fail'))"

# Build {{{1
################################################################

$(ZIPPED_PKG) build: doc
	R CMD build .

# Documentation {{{1
################################################################

doc:
	R -q -e "devtools::document('$(CURDIR)')"

vignettes:
	@echo Build vignettes for already installed package, not from local soures.
	R -q -e "devtools::clean_vignettes('$(CURDIR)')"
	R -q -e "devtools::build_vignettes('$(CURDIR)')"

# Clean {{{1
################################################################################

clean:

# PHONY {{{1
################################################################################

.PHONY: all clean check
