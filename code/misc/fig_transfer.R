# copy over figures for mike

rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
setwd(basedir)
injection.site <- c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl')
injection.site.label <- paste0(injection.site,collapse='-')
c.max = 0.2
opdir <- paste('TauDiffusion032520_Inject',injection.site.label,'_CMax',c.max,'/',sep='')
fig.dir <- paste0(basedir,'figures/')
dir.create(fig.dir,recursive = T)

# Figure 1&2: bidirectional model fits for NTG and G20
f.defs <- list(list(grp='NTG',f='Figure1'),list(grp='G20',f='Figure2'))
for(f.def in f.defs){
  f <- f.def$f
  grp <- f.def$grp
  fig.dir.j <- paste0(fig.dir,f,'/') # make figure directory
  dir.create(fig.dir.j,recursive = T)
  fig.dir.j <- paste0(fig.dir.j,f) # prefix for figure direction/Figure
  desc <- c(grp,'(a) A combination of retrograde and anterograde diffusion models explains pathology spread.',
            '(b) The same as (a), coloring each point by hemisphere relative to injection site.',
            '(c) Independent contributions of anterograde and retrograde spread over time, measured by standardized regression betas. Variance Inflation Factor reveals multicollinearity is not an issue',
            '(d) Vulnerability/residuals from (a) vs MAPT expression.')
  write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'CNDRSpaceFit_bidirectionaladditivemodel.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'CNDRSpaceFitHemiColor_bidirectionaladditivemodel.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'CNDRSpaceFitAnteroRetroBetas_bidirectionaladditivemodel.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'CNDRSpaceVulnerability_bidirectionalmodel_vsMaptRP_071204_01_D02.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'vulnerability_bidirectional_hemiaverage_exclude1 MPI.csv')
  file.copy(from=init,to=paste0(fig.dir.j,'a_hemiaveragevulnerabilityexclude1MPI.csv'))
  init <- paste0(opdir,'diffmodel/bidirectional/',injection.site.label,'/',grp,'log10predictedpath_bidirectional.csv')
  file.copy(from=init,to=paste0(fig.dir.j,'a_log10predicted.csv'))
}

# Figure 3. Retrograde alone: CNDRSpaceCSweepByTimePoint_basemodel.pdf
# shows that average time constant is a good approximation of "real" time constant, 1 MPI not biasing significantly
fig.dir.j <- paste0(fig.dir,'Figure3','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure3') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) For retrograde model, average time constant is a good approximation of "real" time constant, 1 MPI not biasing significantly. vertical line is group average time constant.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/',injection.site.label,'/',grp,'CNDRSpaceCSweepByTimePoint_basemodel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))

# Figure 4. Null model of euclidean distance as network
fig.dir.j <- paste0(fig.dir,'Figure4','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure4') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) Predicted vs. actual for optimized diffusion model where the network is the matrix of inverse euclidean distances between the center of mass of each ABA region.',
          '(b) Same plot as (a), only I have excluded the outliers where predicted path > 95th %ile. Empirically this always includes the injection sites')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean_ExcludeInjectionSites.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))

# in silico injections
file.copy(from=paste0(opdir,'insilico_injections'),to=paste0(fig.dir),recursive = T)

# Figure 5. Relationship between vulnerability and pathology in NTG-G20, and consistency of vulnerability across hemispheres and time.
fig.dir.j <- paste0(fig.dir,'Figure5','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure5') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c('(a) ratio of G20 to NTG regional pathology plotted against time-point specific vulnerability',
          '(b) ratio of G20 to NTG regional pathology plotted against hemisphere and time averaged vulnerability, excluding 1 MPI',
          grp,'(c) Vulnerability values compared between hemispheres at each time point using wilcox rank sum. p-vals not adjusted for MC.',
          '(d) Spatial similarity of model residuals between each hemisphere-time point combo.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'G20vsNTG/path_vs_vuln/bidirectional/',injection.site.label,'/G2019-NTGvsNTGVulnerability_Bidirectional_HemiTimeAverage_exclude1 MPI.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'G20vsNTG/path_vs_vuln/bidirectional/',injection.site.label,'/G2019-NTGvsNTGVulnerability_Bidirectional_TimeDependent.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'diffmodel/vuln_time_hemi/',injection.site.label,'/',grp,'CompareVulnerabilityBetweenHemispheres.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'diffmodel/vuln_time_hemi/',injection.site.label,'/',grp,'VulnerabilityByTimeAndHemisphere.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))

# alternate seed sites
fig.dir.j <- paste0(fig.dir,'Figure6','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure6') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) distance between random seed sites is within 10% of distance between tested injection sites',
          '(b) alternate seeds fit worse than real seeds for 3-9 MPI',
          '(c) Alternate seed fit is partially explained by in projection similarity, out projection similarity, distance to injection site,',
          'and distance between alternate seed sites')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/AlternateSeedSiteDistances.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'SeedSpecificity.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'RandomSeedFitsVsConnectivity.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
