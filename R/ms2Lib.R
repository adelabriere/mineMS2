#' @include LossFormula.R


###CORRESPONDANCE TABLE
CORR_TABLE <- list("D"="dags","S"="spectra","P"="patterns","L"="losses")#"F"="fragments"
UNKNOWN_FORMULA <- NA_character_



###Accessors
setMethod("mm2Spectra","ms2Lib",function(m2l){
	return(m2l@spectra)
})

setMethod("mm2SpectraInfos","ms2Lib",function(m2l){
	return(m2l@spectraInfo)
})

setMethod("mm2Dags","ms2Lib",function(m2l){
	return(m2l@dags)
})

setMethod("mm2Ids","ms2Lib",function(m2l){
	return(m2l@ids)
})

setMethod("mm2Atoms","ms2Lib",function(m2l){
  return(m2l@atoms)
})

setMethod("mm2EdgesLabels","ms2Lib",function(m2l){
	return(m2l@losses)
})

setMethod("mm2NodesLabels","ms2Lib",function(m2l){
	return(m2l@fragments)
})

setMethod("mm2Patterns","ms2Lib",function(m2l){
	return(m2l@patterns)
})

setMethod("mm2Lattice","ms2Lib",function(m2l){
	return(m2l@lattice)
})

setMethod("mm2LatticeIdxItems","ms2Lib",function(m2l){
	return(m2l@latticeIdxItems)
})

setMethod("mm2ReducedPatterns","ms2Lib",function(m2l){
	return(m2l@reducedPatterns)
})

setMethod("mm2ReducedLattice","ms2Lib",function(m2l){
	return(m2l@reducedLattice)
})


###Setter
setMethod("mm2Spectra<-","ms2Lib",function(m2l,value){
	m2l@spectra <- value
	m2l
})

setMethod("mm2SpectraInfos<-","ms2Lib",function(m2l,value){
	m2l@spectraInfo <- value
	m2l
})

setMethod("mm2Ids<-","ms2Lib",function(m2l,value,check=TRUE){
	if(check & any(startsWith(value,names(CORR_TABLE)))){
		stop("Forbidden prefixes for ids: ",paste(names(CORR_TABLE),collapse = ", "),
			 " found in ",value[startsWith(value,names(CORR_TABLE))])
	}

	if(length(value) != length(mm2Spectra(m2l))){
		stop("Number of furnished ids (",paste(length(value)),
			 ") should equal to spectra number (",paste(length(mm2Spectra(m2l))),")")
	}
	m2l@ids <- value
	###Check of the correctness of the IDs.
	m2l
})

#' @export
setMethod("setIds","ms2Lib",function(m2l,ids){
	mm2Ids(m2l) <- ids
	m2l
})



setMethod("mm2Dags<-","ms2Lib",function(m2l,value){
	m2l@dags <- value
	m2l
})

setMethod("mm2EdgesLabels<-","ms2Lib",function(m2l,value){
	m2l@losses <- value
	m2l
})

setMethod("mm2NodesLabels<-","ms2Lib",function(m2l,value){
	m2l@fragments <- value
	m2l
})

setMethod("mm2Patterns<-","ms2Lib",function(m2l,value){
	m2l@patterns<- value
	m2l
})

setMethod("mm2Lattice<-","ms2Lib",function(m2l,value){
	m2l@lattice<- value
	m2l
})

setMethod("mm2LatticeIdxItems<-","ms2Lib",function(m2l,value){
	m2l@latticeIdxItems<- value
	m2l
})

setMethod("mm2ReducedPatterns<-","ms2Lib",function(m2l,value){
	m2l@reducedPatterns<- value
	m2l
})

setMethod("mm2ReducedLattice<-","ms2Lib",function(m2l,value){
	m2l@reducedLattice<- value
	m2l
})

setMethod("mm2Atoms<-","ms2Lib",function(m2l,value){
  m2l@atoms<- value
  m2l
})


isLoss <- function(m2l){
	m2l@loss
}

###A list of the recognised format
#' @export
recognisedFormat <- function(){
	return(c("mgf"))
}


###Parse an mgf and return a list of the Spectrum2 object
parse_mgf_spectrum2 <- function(filename){
	msnexp <- MSnbase::readMgfData(filename)
	lspec <- vector(mode="list",length=length(msnexp))
	for(i in 1:length(msnexp)){

		lspec[[i]] <- msnexp[[paste("X",i,sep="")]]
	}
	
	metadata <- fData(msnexp)[paste("X",1:length(msnexp),sep=""),]
	rlist <- list(spec=lspec,supp=metadata)

	return(rlist)
}

###Check the format fo a list of files.
checkFormat <- function(x){
	rf <- recognisedFormat()
	splitted <- strsplit(x,".",fixed=TRUE)

	exts <- sapply(splitted,tail,1)

	is_ok <- exts %in% rf

	if(any(!(is_ok))) stop(paste("Unrecognized format:",exts[!is_ok],"in",x[!is_ok]))

	return(exts)
}

parseMS2file_line <- function(x){
	do.call(paste("parse",x[2],"spectrum2",sep="_"),list(filename=x[1]))
}

make_initial_title <- function(spec_infos){
	mz_str <- sprintf("%0.3f",spec_infos[,"mz.precursor"])
	paste("Precursor : ",mz_str," (S",1:nrow(spec_infos),")",sep="")
}


convert_formula <- function(form_vec){
  ###We remoe the whitespaces and the special character
  form_vec <- trimws(str_replace(form_vec,"\\?|\\-|\\+|\\-",""))
  
  pnf <- which(is.na(form_vec)|(nchar(form_vec)==0))
  form_vec[pnf] <- UNKNOWN_FORMULA
  
  
  form_vec
}



#' ms2Lib constructor
#'
#' Create a ms2Lib object from file sin the correct format. Eventually add some supplementary informations given in a file. These supplementary informations
#' may notably include the composition or the ofrmula of the molecules, which will then be used to plot the created graphs..
#'
#' @export
#' @param x May be one of the following
#' \itemize{
#' \item A character vector giving the path to a drictory full of readable format.
#' \item A list of spectrum2 object which will be integrated directly.
#' \item A single .mgfspectrum regrouping multiple files.
#' }
#' @param suppInfos Supplementary information to be associated to the spectra.
#' It should be of the same size as the number of spectra. If there is a "file" column, this
#' column is used to match the filenames. A "composition" or "formula" fields may also be present, and is then used to set the
#' formula of the molecules.
#' @param ids A supplementary vector giving a set of ids to design the spectra. It may be any character vector which
#' does not start with \code{'P', 'L', 'S'} as they are used internally by mineMS2. Alternatively if a suppInfos table is furnished
#' and it contains an id fields, it will be used. If no ids are furnished, an id will be generated for each spectra in the form \code{'S1', 'S2', ..., 'SN'}
#' where N is the number of furnished spectra.
#' @param infosFromFiles Shall the other informations present in the files be added to the supplementary infos.
#' @examples
#' print("examples to be put here")
ms2Lib <- function(x, suppInfos = NULL,ids = NULL, intThreshold = NULL, infosFromFiles = FALSE){

	m2l <- new("ms2Lib")

	origin <- "R"
	lfiles <- NULL

	suppMetadata <- NULL
	###The kind of the acquisition is assessed there.
	if(class(x)=="list"){
		if(all(sapply(x,class) == "Spectrum2")){
			mm2Spectra(m2l) <- x
		}else{
			message("Unrecognized input, use one of the constructor described in the ms2Lib doc.")
		}

	}else if(class(x)=="character"){
		origin <- "file"
		if(length(x)==1){
		###Single mgf file or dir
			if(dir.exists(x)){
				lfiles <- list.files(x,full.names = TRUE)
				exts <- checkFormat(lfiles)
				message("Reading ",length(exts)," files with format(s): ",unique(exts))
				
				tres <- apply(matrix(c(lfiles,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)
				mm2Spectra(m2l) <- sapply(tres,"[",i="spec")
				suppMetadata <- sapply(tres,"[",i="supp")
				
			}else{ ###Case of a single spectra.
				message("Reading MGF file ", x, ".")
				tres <- parseMS2file_line(c(x, 'mgf'))
				mm2Spectra(m2l) <- tres$spec
				suppMetadata <- tres$supp
				
			}
		}else{
			###Case of multiples singles spectra
			exts <- checkFormat(x)
			message("Reading ",length(exts)," files with format(s): ",unique(exts))
			mm2Spectra(m2l) <- apply(matrix(c(x,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)
		}
	}

	mm2Spectra(m2l) <- do.call("c",mm2Spectra(m2l))
	
	if(length(mm2Spectra(m2l))>1000){
	  stop("At the moment it is impossible ot process more than 1000 spectra at the same time.")
	}


	###data.frame is initialized. With the mass of the precusors
	temp_df <- data.frame("mz.precursor" = sapply(mm2Spectra(m2l),function(x){
												  precursorMz(x)}))


	if(!is.null(intThreshold)){
		message("Removing peaks with an intensity lower than ",intThreshold)
		for(is in seq_along(mm2Spectra(m2l))){
			m2l@spectra[[is]] <- removePeaks(m2l@spectra[[is]],t = intThreshold)
			m2l@spectra[[is]] <- clean(m2l@spectra[[is]],all=TRUE)
		}
	}

	if(origin=="file"){
		if(is.null(lfiles)){
			if(length(x)!=nrow(temp_df)){
				temp_df$file <- rep(x,nrow(temp_df))
			}else{
				temp_df$file <- x
			}
		}else{
			temp_df$file <- lfiles
		}
	}

	temp_df$title <- make_initial_title(temp_df)

	##Adding the supplementary information if necessary while check for an id fileds.
	if(!is.null(suppInfos)){
		if(nrow(suppInfos)!= length(m2l@spectra)){
			stop("The number of suppInfos rows (",nrow(suppInfos),
				 ") do not match the number of spectra (",
				 length(m2l@spectra),")furnished")
		}else{
			if("file" %in% colnames(suppInfos)){
				pm <- match(temp_df$file,suppInfos$file)
				if(any(is.na(pm))) stop('"file" column furnished, but there was an error matching it against files.')

				temp_df <- cbind(temp_df,suppInfos[pm,])
			}else{
				temp_df <- cbind(temp_df,suppInfos)
			}

			if(("id" %in% colnames(suppInfos)) &
			   (is.null(ids))){
				mm2Ids(m2l) <- suppInfos[,"id"]
			}
		}

	}else{
		if(!is.null(ids)){
			mm2Ids(m2l) <- ids
		}else{
			mm2Ids(m2l,check=FALSE) <- paste("S",1:length(mm2Spectra(m2l)),sep="")
		}
	}
	
	if(infosFromFiles&!is.null(suppMetadata)){
	  temp_df <- cbind(temp_df,suppMetadata)
	}
	
	
	####Adding the molecular formula.
	cnames <- tolower(colnames(temp_df))
	
	pf <- which(cnames %in% c("formula","composition"))
	if(length(pf)==0){
	  message("No 'formula' column found. All formula are considered as unknown.")
	  temp_df$formula <- rep(UNKNOWN_FORMULA,nrow(temp_df))
	}else{
	  temp_df[,pf] <- convert_formula(as.character(temp_df[,pf]))
	}

	mm2SpectraInfos(m2l) <- temp_df
	m2l
}


get_formula <- function(m2l){
  vf <- match("formula",tolower(colnames(m2l@spectraInfo)))
  return(m2l@spectraInfo[,vf])
}

#' @export
setMethod("show","ms2Lib",function(object){
	cat("An ms2Lib object containing",length(object),"spectra.\n")
  if(length(mm2Atoms(object))!=0){
	  cat("It has",nrow(mm2EdgesLabels(object)),"edges labels built with atoms",paste(mm2Atoms(object),collapse=","),".\n")
  }
  cat("The available supplementary informations are:",colnames(mm2SpectraInfos(object)),"\n")
	cat("It contains: ",length(mm2Patterns(object)),"patterns\n")
	if(length(mm2ReducedPatterns(object))!=0) cat("It has been reduced to ",
											   length(mm2ReducedPatterns(object)),"patterns")
})

#' @export
setMethod("length","ms2Lib",function(x){
	return(length(mm2Spectra(x)))
})

getFormula <- function(m2l){
  vnames <- tolower(colnames(m2l@spectraInfo))
  return(as.character(m2l@spectraInfo[,match("formula",vnames)]))
}



#' Mine recurrent subgraph from a set of graphs.
#'
#' Mine all the complete recurring subgraphs.
#'
#' @param m2l An m2Lib object to be processed.
#' @param count The number of spectra in which the spectrum need to be sampled.
#' @param sizeMin The minimum size of the mined patterns.
#' @param kTree The maximum depth of the path tree.
#' @param precursor Should only the occurences coming form the root be conserved.
#'
#' @return The filled ms2Lib object.
#' @export
#'
#' @examples
#' print("examples to be put here")
setMethod("mineClosedSubgraphs","ms2Lib",function(m2l, count = 2, sizeMin = 2,kTree = NULL, precursor = FALSE){
	if(count<2){
		warning("'count' parameters set to ",count," it is therefore set to 2.")
		count <- 2
	}


	###Get the data.frame correspoding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No labels constructed, use the DiscretizeMallLosses function first.")
	}

	if(is.null(kTree)){
	if(nrow(mm2EdgesLabels(m2l))<600){
		kTree <- 3
	}else{
		kTree <- 2
	}
	}



	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("sizeMin parameters set to ",sizeMin," risk of computational overhead.")
	}

	###Converting the dags into data.frame.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)

	###Mining the patterns.
	resRcpp <- mineClosedDags(df_vertices,df_edges,processing,count,kTree,sizeMin,precursor)

	mm2LatticeIdxItems(m2l) <- resRcpp$items

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- sapply(resRcpp$patterns,fragPattern,USE.NAMES = FALSE)

	###Initializing the names of the patterns.
	for(i in 1:length(m2l@patterns)) mm2Name(m2l@patterns[[i]]) <- paste("P",i,sep="")


	###Buillding of the lattice
	mm2Lattice(m2l) <- graph_from_data_frame(resRcpp$edges,directed=TRUE,resRcpp$nodes)

	###We add a label filed which give all the values of this.


	message("Processing finished, ",length(mm2Patterns(m2l))," patterns mined.")
	m2l
})



setMethod("mineClosedSubgraphs","ms2Lib",function(m2l, count = 2, sizeMin = 2, precursor = FALSE){
	if(count<2){
		warning("'count' parameters set to ",count," it is therefore set to 2.")
		count <- 2
	}


	###Get the data.frame correspoding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No labels constructed, use the DiscretizeMallLosses function first.")
	}


	kTree <- NULL
	if(nrow(mm2EdgesLabels(m2l))<600){
		kTree <- 2
	}else{
		kTree <- 1
	}



	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("sizeMin parameters set to ",sizeMin," risk of computational overhead.")
	}

	###Converting the dags into data.frame.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)

	###Mining the patterns.
	resRcpp <- mineClosedDags(df_vertices,df_edges,processing,count,kTree,sizeMin,precursor)

	mm2LatticeIdxItems(m2l) <- resRcpp$items

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- sapply(resRcpp$patterns,fragPattern,USE.NAMES = FALSE)

	###Initializing the names of the patterns.
	for(i in 1:length(m2l@patterns)) mm2Name(m2l@patterns[[i]]) <- paste("P",i,sep="")


	###Buillding of the lattice
	mm2Lattice(m2l) <- graph_from_data_frame(resRcpp$edges,directed=TRUE,resRcpp$nodes)

	###We add a label filed which give all the values of this.


	message("Processing finished, ",length(mm2Patterns(m2l))," patterns mined.")
	m2l
})



###Parse an id
parseId <- function(m2l,idx){

	prefix <- substring(idx,1,1)
	number <- as.integer(substring(idx,2))
	if(!(prefix %in% names(CORR_TABLE))){
		###Checking if it's in the ids fields.
		resm <- match(idx,mm2Ids(m2l))
		if(is.na(resm)){
			stop("Invalid prefix ",prefix," authorized prefix are ",
			 	paste(names(CORR_TABLE),collapse=", "))
		}else{
			return(list(type=CORR_TABLE[["S"]],num=resm))
		}
	}

	###The case of the L prfix is handled directly.
	if(prefix=="L"){
		if(nrow(mm2EdgesLabels(m2l))<number) stop("Invalid id for mass_losses: ",number,".")
		return(list(type=CORR_TABLE[[prefix]],num=number))
	}
	if( (number<=length(slot(m2l,CORR_TABLE[[prefix]])))&
		(number>=1)){
		return(list(type=CORR_TABLE[[prefix]],num=number))
	}else{
		stop("Invalid id for ",CORR_TABLE[[prefix]]," : ",number,".")
	}
}

mm2get <- function(m2l,arglist){
	if(class(slot(m2l,arglist[[1]]))=="data.frame"){
		return(slot(m2l,arglist[[1]])[arglist[[2]],])
	}
	(slot(m2l,arglist[[1]]))[[arglist[[2]]]]
}

#' Return the dag correspodning to a spectra or the modif.
#'
#' Indexing function.
#'
#' @param x An ms2Lib oject.
#' @param i The index, a string starting by S if it a spectrum, P if it's a pattern D if it's a dag, L if it's a loss
#' F if it's a fragment.
#' @param drop
#'
#' @return A list containg the object.
#' @export
#'
#' @examples
setMethod('[','ms2Lib',function(x,i,j=NULL,...,drop=TRUE){
	if(length(i)==1){
		 temp <- mm2get(x,parseId(x,i))
		 # if(drop){
		 	return(temp)
		 # }
	}else{
		res <- lapply(lapply(i,parseId,m2l=x),mm2get,m2l=x)
		return(res)
	}
})


#' Plot an element given an idx
#'
#' The method depends of the type of the ID furnished. The following prefixes are supported:
#' \code{'P'} A pattern is called plot method of fragPattern object.
#' \code{'S'} A spectrum is plotted calling the Plot method of spectrum 2 object.
#' Any other value will be removed.
#'
#' @param x An ms2Lib oject.
#' @param y The index, a string starting by S if it a spectrum, P if it's a pattern or D if it's a dag.
#' @param ... supplementary arguments to be passed by the method.
#'
#' @return a fragPattern object of an igraph graph object.
#' @export
#'
#' @examples
setMethod("plot", "ms2Lib",
		  function(x,
		  		 y,title=NULL,
		  		 ...) {
		  	if(length(y)>1){
		  		warning("A single if may be plotted on each call, plotting the first element only")
		  		y <- y[1]
		  	}
		  	rid <- parseId(x,y)
		  	if(rid[[1]]=="patterns"){
		  	  toccs <- x[y]@occurences[,1]
				return(plot(x[y],title = y,dags=mm2Dags(x),edgeLabels=(mm2EdgesLabels(x)),
				     atoms=names(x@atoms),formula=get_formula(x)[toccs],...))
		  	}else if(rid[[1]]=="spectra"){
		  		plot_Spectrum2(x[y],full=TRUE,...)
		  	}else if(rid[[1]]=="dags"){
		  	  if(is.null(title)) title="Fragmentation Graph"
		  		plot_dag(x[y],idx=y,edgeLabels=(mm2EdgesLabels(x)),atoms=x@atoms,title=title,...)
		  		# stop("DAGS plotting not implemented at the moment.")
		  	}else if(rid[[1]]=="losses"){
		  		stop("Impossible to plot a loss")
		  	}else if(rid[[1]]=="fragments"){
		  		stop("Impossible to plot a fragment")
		  	}
})

#
# setMethod('[[','MSMSacquisition',function(x,i,j,...,drop=TRUE){


####Search and info functions
findMz.S <- function(m2l,mz,tol){
	infos <- mm2SpectraInfos(m2l)
	matched <- which(abs(infos$mz.precursor-mz)<tol)
	if(length(matched)==0) return(character(0))
	paste("S",matched,sep="")
}

findMz.L <- function(m2l,mz,tol){
	infos <- mm2EdgesLabels(m2l)
	matched <- which(abs(infos$mz-mz)<tol)
	if(length(matched)==0) return(character(0))
	paste("L",matched,sep="")
}

#' Search in an ms2Lib object?
#'
#' Search an ms2Lib object givena tolerance in ppm or dmz.
#'
#' @param m2l ms2Lib object
#' @param mz A double giving the mass to be searched.
#' @param type The data to be search, "S" for spectra and "L" for losses
#' @param ppm The tolerance in ppm
#' @param dmz The minimum tolerance in Da, if the ppm tolerance is lower in Da than this threshold, this threshold is selected.
#'
#' @return A character vector giving the IDs of the found losses or elements.
#' @export
#'
#' @examples
#' print("examples to be put here")
findMz <- function(m2l,mz,type=c("S","L"),ppm=15,dmz=0.01){
	if(class(m2l)!="ms2Lib") stop("m2l should be an 'ms2Lib' object.")
	type <- match.arg(type)
	tol <- max(dmz,mz*ppm*1e-6)
	if(type=="S"){
		if(length(mz)>1){
			return(sapply(mz,findMz.S,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.S(m2l,mz,tol))
	}else{
		if(length(mz)>1){
			return(sapply(mz,findMz.L,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.L(m2l,mz,tol))
	}
}


getInfo.L <- function(num,m2l){
	titles <- colnames(mm2EdgesLabels(m2l))
	titles <- titles[!(titles %in% 	c("sig", "fused",
									  "adv_loss", "pen_loss", "carb_only", "nitrogen_only", "ch_only",
									  "full_labels", "labs"))]




	return(mm2EdgesLabels(m2l)[
		num,titles,drop=FALSE])

}

getInfo.S <- function(num,m2l){
	titles <- colnames(mm2SpectraInfos(m2l))
	titles <- titles[!(titles %in% 	c("title"))]
	return(mm2SpectraInfos(m2l)[
		num,c(titles),drop=FALSE])
}

#' Return the available informationon a lost of a spectra.
#'
#' @param m2l An ms2Lib object
#' @param ids A vector of IDs should be a spectrum or a losses.
#' @param ... Supplementary information ot be passed to the getInfo function.
#'
#' @return A data.frame giving informations about the queried elements.
#' @export
#'
#' @examples
#' print("Examples to be put here")
getInfo <- function(m2l,ids){
	if(class(m2l)!="ms2Lib") stop("m2l should be an 'ms2Lib' object.")
	authorizedValue <- c("losses","spectra")

	if((length(ids)==1) && (nchar(ids)==1)){
		if(ids=="S"){
			ids <- vrange(m2l,"S")
		}else if(ids=="L"){
			ids <- vrange(m2l,"L")
		}else{
			stop("Single letters are only allowed for S and L.")
		}
	}

	pids <- sapply(ids,parseId,m2l=m2l,simplify=FALSE)

	type <- sapply(pids,'[[',i="type")

	if(all(type %in% authorizedValue)){
		num <- sapply(pids,'[[',i=2)
		res <- sapply(pids,function(x,m2l){
			if(x[[1]] == "losses"){
				return(getInfo.L(x[[2]],m2l))
			}
			if(x[[1]] == "spectra"){
				return(getInfo.S(x[[2]],m2l))
			}
		},simplify=FALSE,m2l=m2l)
		if(length(unique(type))==1){
			return(do.call("rbind",res))
		}else{
			return(res)
		}

	}else{
		stop("Invalid type for getInfo: ",unique(type[!(type %in% authorizedValue)]))
	}
}


#' Return the range of iteration on an ms2Lib object.
#'
#' @param m2l AN ms2Lib object
#' @param type "S","L","P" or "F"
#' @param reduced Used only if "type" is set to "P".
#'
#' @return A character vector giving the existing ids.
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("vrange","ms2Lib",function(m2l,type=c("S","L","P","F"), reduced=TRUE, as.number=FALSE){
	type <- match.arg(type)
	if(type=="S"){
		if(length(m2l)==0) return(character(0))
	  if(as.number) return(1:length(m2l))
		return(paste("S",1:length(m2l),sep=""))
	}
	if(type=="P"){
		if(length(mm2Patterns(m2l))==0) return(character(0))
		if(reduced){
		  if(as.number) return(mm2ReducedPatterns(m2l))
			return(paste("P",mm2ReducedPatterns(m2l),sep=""))
		}else{
		  if(as.number) return(1:length(mm2Patterns(m2l)))
			return(paste("P",1:length(mm2Patterns(m2l)),sep=""))
		}
	}
	if(type=="L"){
		if(nrow(mm2EdgesLabels(m2l))==0) return(character(0))
		return(paste("L",1:nrow(mm2EdgesLabels(m2l)),sep=""))
	}
	if(type=="F"){
		if(nrow(mm2NodesLabels(m2l))==0) return(character(0))
	  if(as.number) return(1:nrow(mm2NodesLabels(m2l)))
		return(paste("L",1:nrow(mm2NodesLabels(m2l)),sep=""))
	}
})


hasPatterns <- function(m2l){
  return(length(m2l@patterns)!=0)
}

setMethod("hasCoverage","ms2Lib",function(x){
  if(hasPatterns(x)){
    return(COVERAGE_NAME %in% colnames(x@patterns[[1]]@occurences))
  }else{
    stop(paste("Patterns need to be computed before obtaining coverage.",sep=""))
  }
  return(FALSE)
})


#' Calculate the coverage of all the patterns in the dataset.
#'
#' @param m2l The ms2Lib to bo computed.
#'
#' @return The m2l object with all the coverage calculated.
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("calculateCoverage","ms2Lib",function(x){
  loss_mz <- mm2EdgesLabels(x)$mz
  mgs <- mm2Dags(x)
  pb <- txtProgressBar(min = 0, max = length(mm2Patterns(x)), initial = 0, char = "=",
                 width = NA, "Covergae calculation", "cov_calc", style = 3, file = "")
  for(i in seq_along(mm2Patterns(x))){
    setTxtProgressBar(pb, i)
    x@patterns[[i]] <- calculateCoverage(x@patterns[[i]],loss_mz,mgs)
  }
  return(x)
})
