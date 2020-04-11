source('code/misc/miscfxns.R')

lm.mask.ant.ret <- function(y,Xt.retro,Xt.antero){
  # fit linear model: y = b0 + b.r*Xt.retro + b.a*Xt.antero + error
  # given Xt.antero and Xt.retro for a given set of time constants
  # return pearson r between y and y-hat (predicted vs actual path)
  #
  # INPUTS:
  # y: observed log pathology
  # Xt.retro: predicted pathology from retrograde model
  # Xt.antero: predicted pathology from retrograde model
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  Xt.retro <- log(Xt.retro,base=10)
  Xt.antero <- log(Xt.antero,base=10)
  mask <- y == -Inf | Xt.retro == -Inf | Xt.antero == -Inf | is.na(Xt.retro) | is.na(Xt.antero)
  df.m <- data.frame(y=y[!mask],x1=Xt.retro[!mask],x2=Xt.antero[!mask])
  m <- lm(y~x1+x2,data=df.m)

  return(cor(m$fitted.values,df.m$y))
}

c.CNDRspace.objective <- function(params.opt,log.path,tps,L.out.retro,L.out.antero,Xo,ABA.to.CNDR.key,fxn){
  # fits time constant by modeling CNDR data
  # INPUTS:
  # OPTIMIZED:
  # params.opt: parameters to optimize cor(observed path, b0 + b.r*retro(c.retro) + b.a*antero(c.antero))
  #     -c.retro: time constant for retrograde model --- params[1]
  #     -c.antero: time constant for anterograde model --- params[2]
  #     -lm.mask.ant.ret fits linear model to find b0, b.a, b.r for each c.retro,c.antero pair
  #     -this is done separately from optimization and params are discarded
  #
  # STATIC:
  # log.path: list of vectors of log-10 transformed pathology scores *in CNDR space*  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out.retro: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads axon to dendrite (retrograde)
  # L.out.antero: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads dendrite to axon (anterograde)
  # Xo: vector of initial pathology
  # ABA.to.CNDR.key: key to convert ABA to CNDR regions
  # fxn: matrix exponential function (fastest is scipy.linalg through reticulate)
  
  # this function runs diffusion model in ABA space, 
  # but computes correlation with real data in CNDR annotation space to assess fit of time constant c
  # CNDR names are names of log.path, ABA names are names of Xo
  
  # OUTPUTS:
  # r: mean pearson r value across all provided time points for given set of params
  
  ####
  
  # distribute parameters to interpretable names 
  c.retro <- params.opt[1]
  c.antero <- params.opt[2]
  
  ptm <- proc.time()
  Xt.retro <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.retro,Xo,c.retro,t.,fxn=scipy.linalg$expm))) # predict path using linear diff model into matrix that is region-by-time
  Xt.antero <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.antero,Xo,c.antero,t.,fxn=scipy.linalg$expm))) # predict path using linear diff model into matrix that is region-by-time
  print(paste0('predict: ',(proc.time()-ptm)['elapsed']))
  print(paste0('optim,',c.retro,',',c.antero)) # output tested parameters to see what optim is doing
  #ptm <- proc.time()
  Xt.retro <- quiet(map.ABA.to.CNDR(Xt.retro,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  Xt.antero <- quiet(map.ABA.to.CNDR(Xt.antero,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  #print(paste0('map: ',(proc.time()-ptm)['elapsed']))
  Xt.sweep <- sapply(1:length(tps), function(t.) lm.mask.ant.ret(log.path[[t.]],Xt.retro[,t.],Xt.antero[,t.])) # this replaces the two lines below
  
  print(mean(Xt.sweep)) # display value of fit/objective function
  return(mean(Xt.sweep)) # return mean across provided tps
}


DEPRECATED.c.CNDRspace.objective <- function(params,log.path,tps,L.out.retro,L.out.antero,Xo,ABA.to.CNDR.key){
  # DEPRECATED: realized better method partway through but saving temporarily in case i need to return to this
  # fits time constant by modeling CNDR data
  # INPUTS:
  # OPTIMIZED:
  # params: parameters to optimize cor(observed path, b0 + b.r*retro(c.retro) + b.a*antero(c.antero))
  #     -c.retro: time constant for retrograde model --- params[1]
  #     -c.antero: time constant for anterograde model --- params[2]
  #     -b0: intercept --- params[3]
  #     -b.r: beta to weight retrograde prediction --- params[4]
  #     -b.a: beta to weight anterograde prediction --- params[5]
  #     **** another thing to do is: cor(observed path, b0 + b.r*retro(c.retro) + (1-b.r)*antero(c.antero))
  #
  # STATIC:
  # log.path: list of vectors of log-10 transformed pathology scores *in CNDR space*  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out.retro: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads axon to dendrite (retrograde)
  # L.out.antero: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads dendrite to axon (anterograde)
  # Xo: vector of initial pathology
  # ABA.to.CNDR.key: key to convert ABA to CNDR regions
  
  # this function runs diffusion model in ABA space, 
  # but computes correlation with real data in CNDR annotation space to assess fit of time constant c
  # CNDR names are names of log.path, ABA names are names of Xo
  
  # OUTPUTS:
  # r: mean pearson r value across all provided time points for given set of params
  
  ####
  
  # distribute parameters to interpretable names 
  c.retro <- params[1]
  c.antero <- params[2]
  b0 <- params[3]
  b.r <- params[4]
  b.a <- params[5]
  
  scipy.linalg <- reticulate::import('scipy.linalg') # import python function for matrix exponential because it's faster
  #log.path <- lapply(1:length(log.path), function(t) log.path[[t]][mask[[t]]])
  # compute fit at each time point for range of time
  Xt.sweep <- matrix(NA,nrow=1,ncol=length(tps))
  print(paste0('c ',c.i,' out of ',length(c.rng)))
  ptm <- proc.time()
  Xt.retro <- do.call('cbind',lapply(tps,function(t) predict.Lout(L.out.retro,Xo,c.val,t,fxn=scipy.linalg$expm))) # predict path using linear diff model into matrix that is region-by-time
  print(paste0('predict: ',(proc.time()-ptm)['elapsed']))
  #ptm <- proc.time()
  Xt.c <- quiet(map.ABA.to.CNDR(Xt.c,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  #print(paste0('map: ',(proc.time()-ptm)['elapsed']))
  Xt.sweep[c.i,] <- sapply(1:length(tps), function(t) cor.mask(Xt.c[,t],log.path[[t]])) # this replaces the two lines below
  
  return(Xt.sweep)
}