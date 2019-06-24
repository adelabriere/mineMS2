# vi: fdm=marker ts=4 et cc=80 tw=80

# Constants {{{1
################################################################

TEST.DIR <- file.path(getwd(), '..')
RES.DIR  <- file.path(TEST.DIR, 'res')

# Test basic processing {{{1
################################################################################

test_basic_processing <- function() {
    
    path_mgf <- file.path(RES.DIR, 'ex_mgf.mgf')

    m2l <- mineMS2::ms2Lib(path_mgf)
    infos <- mineMS2::getInfo(m2l, "S")
    ids <- paste(paste0("MZ", infos[,"mz.precursor"]), sep = "_")
    m2l <- mineMS2::setIds(m2l, ids)
    m2l <- mineMS2::discretizeMassLosses(m2l,
                                         dmz = 0.008,
                                         ppm = 8,
                                         heteroAtoms = FALSE,
                                         maxFrags = 15)
    m2l <- mineMS2::mineClosedSubgraphs(m2l,
                                        sizeMin = 1,
                                        count = 2)
    mineMS2::plotPatterns(m2l, full = TRUE)
}

# Main {{{1
################################################################################

testthat::context('Main tests.')
testthat::test_that('A basic processing works fine.', test_basic_processing())
