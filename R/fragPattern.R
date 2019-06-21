#' @include references.R

setMethod("mm2Graph", "fragPattern",function(pat){
	return(pat@graph)
})

setMethod("mm2Occurences", "fragPattern",function(pat){
	return(pat@occurences)
})

setMethod("mm2Root", "fragPattern",function(pat){
	return(pat@root)
})

setMethod("mm2Name", "fragPattern",function(pat){
	return(pat@name)
})


setMethod("mm2Graph<-", "fragPattern", function(pat,value){
	pat@graph <- value
	pat
})

setMethod("mm2Occurences<-", "fragPattern", function(pat,value){
	pat@occurences <- value
	pat
})
setMethod("mm2Root<-","fragPattern", function(pat,value){
	pat@root <- value
	pat
})

setMethod("mm2Name<-", "fragPattern",function(pat,value){
	pat@name <- value
	pat
})


#' @export
setMethod("show","fragPattern",function(object){
	cat("A fragPattern object with",vcount(mm2Graph(object))-1,"losses occuring ",
		nrow(mm2Occurences(object)),"times.\n")
})


#' fragPattern constructor
#'
#' This is an internal constructor which is used to pass form the rcpp
#' output to a fragPattern object including the occurences. It should not be
#' called directly by the user.
#'
#' @param rlist rlist is a list containg 2 fields,
#' \code{edges}: An edge list which includes the "from","to" and "lab" fileds, which by processed
#' by the graph_from_data_frame function.
#' \code{occurences}: A matrix with 2 column, the first giving the id in which the patterns are found
#' and the second giving the nodes.
#'
#' @return
#'
#' @examples
#' print("Exmaples to be put here")
fragPattern <- function(rlist){
	fp <- new("fragPattern")
	mm2Graph(fp) <- graph_from_data_frame(rlist$edges)
	mm2Occurences(fp) <- rlist$occurences

	##The root is always 0
	mm2Root(fp) <- as.integer(0)
	fp
}


setMethod("calculateCoverage","fragPattern",function(x,mzloss,mgs){
  occs <- x@occurences
  coverages <- rep(NA_real_,nrow(occs))
  for(j in 1:nrow(occs)){
    gid <- occs[j,"gid"]
    cg <-mgs[[gid]]
    mapv <- get_mapping(cg ,x@graph, mzloss,root=occs[j,"idx"])
    
    ###Plotting of the spectra
    intv <- vertex_attr(cg, "rel_int")
    mzs <- vertex_attr(cg, "mz")
    ids <- V(cg)
    
    ###Peaks are split between matched and non matched.
    matched_peaks_idx <- match(mapv[2,], ids)
    coverages[j] <- sum(intv[matched_peaks_idx])/sum(intv)
  }
  onames <- colnames(occs)
  if(!COVERAGE_NAME %in% colnames(occs)){
    occs <- cbind(occs,coverages)
    colnames(occs) <- c(onames,COVERAGE_NAME)
  }else{
    occs[,COVERAGE_NAME] <- coverages
  }
  x@occurences <- occs
  x
})

setMethod("hasCoverage","fragPattern",function(x){
  COVERAGE_NAME %in% colnames(x@occurences)
})


relabelOccurrences <- function(pat,newlabs){
	tempOccs <- mm2Occurences(pat)
	tempOccs[,1] <- newlabs[tempOccs[,1]]
	mm2Occurences(pat) <- tempOccs
	pat
}



