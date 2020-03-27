unit.test <- function(test,pass,fail){
	# test: logical expression to evaluate
	# pass: error to print if true
	# fail: error to print if false
	if(test){
		print(pass)
	} else{print(fail)}

}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}

colSD <- function(x){
  return(apply(x,2,sd))
}

fda <- function(x1,x2,nperms = 1000){
	# x1: N-by-P matrix where there are N observational units and P are features, probably time
	# x2: M-by-P matrix where there are N observational units and P are features, probably time
	# perform functional data analysis, where you compute the mean difference between column
	# means of x1 and x2 and compare them to permuted versions of x1 and x2 with the observational units
	# switched	


}

inf.nan.mask <- function(x){
	# INPUTS:
	# x: matrix, df, or vector
	#
	# OUTPUTS:
	# x.masked: x with all rows (if matrix or df) or elements (if vector) containing Infs or NaNs removed
	
	if(is.vector(x)){
		mask <- x %in% c(-Inf,Inf,NaN)# find infs or nans
		x.masked <- x[!mask]
	} else if(is.matrix(x) | is.data.frame(x)){
		mask <- do.call('cbind',lapply(1:ncol(x),function(j) x[,j] %in% c(-Inf,Inf,NaN))) # find infs or nans
		x.masked <- x[rowSums(mask)==0,]
	}
	return(x.masked)
	
}
source.save <- function(script,output){
	# wrapper for source function that saves output and input of script
	file.create(output)
	con <- file(output)
	#sink(con, append=TRUE)
	sink(con, type="output")
	source(script)
	sink() 
	sink(type="output")

}

quiet <- function(x) { 
	# https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message/34208658#34208658
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

name <- function(x,x.names){
	# INPUTS:
	# x: vector or dataframe
	# x.names: names for elements or columns of x
	# OUTPUTS:
	# x with x.names as names
	names(x) <- x.names
	return(x)
}

collapse.columns <- function(df,cnames=colnames(df),groupby=NULL){
  # INPUTS:
  # df: dataframe
  # cnames: column names to perform operation on, default to all columns
  # groupby: column name to group variables by, treated separately from variables in cnames
  
  # OUTPUTS:
  # df.new: dataframe with 2 columns:
  # values: all columns in cnames vertically concatenated. 
  # names: elements of cnames corresponding to rows
  # group: groups of observations in df for each variable in cnames
  
  df.names <- do.call('cbind',lapply(cnames, function(n) rep(n,nrow(as.matrix(df[,cnames])))))  
  df.new <- data.frame(values = as.vector(as.matrix(df[,cnames])),names=as.vector(df.names))
  if(!is.null(groupby)){
    df.grp <- do.call('cbind',lapply(cnames,function(n) df[,groupby]))
    df.new$group <- as.vector(df.grp)
  }
  return(df.new)
}

col.Which.Max <- function(x){
  cwm <- unlist(apply(x,2,function(y) which.max(y)))
  return(cwm)
}

row.Which.Max <- function(x){
  rwm <- unlist(apply(x,1,function(y) which.max(y)))
  return(rwm)
}
