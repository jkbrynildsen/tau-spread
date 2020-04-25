scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default)

params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, lower=c(10e-7,10e-7), # optimization. c's must be > 0
      log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
      Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm) # static inputs