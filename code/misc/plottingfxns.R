plot.Xt <- function(Xt,t){
	# Xt: NxT matrix of nodal value (N) over time (T)
	t <- do.call('rbind',lapply(1:length(Xo), function(i) matrix(t,nrow=1)))
	ROI <- do.call('cbind',lapply(t, function(i) 1:length(Xo)))
	df <- data.frame(y=as.vector(Xt),x=as.vector(t),as.vector(ROI))
	p <- ggplot(df) + geom_line(aes(x=x,y=y,color=ROI))
	p

}

p.xy <- function(x,y,xlab,ylab,ttl='',col='black',alpha=1){
  # INPUTS:
  # x: x variable, vector
  # y: y variable, vector
  # xlab, ylab, ttl: character labels for axes
  # col: point color and line color
  #
  # OUTPUTS:
  # scatter plot with r and p value for pearson correlation between x and y
  # and linear fit

  df <- data.frame(x=x,y=y)
  df <- inf.nan.mask(df)
  c.test <- cor.test(df$x,df$y)
  r.text <- paste0('r = ',signif(c.test$estimate,2),'\np = ',signif(c.test$p.value,2)) # annotation

  p <- ggplot(df) + geom_point(aes(x=x,y=y),color=col,stroke=0,alpha=alpha) + geom_smooth(aes(x=x,y=y),fill=col,color=col,method='lm') + 
    xlab(xlab) + ylab(ylab) + ggtitle(ttl) + 
    annotate("text",size = 2, x = Inf,y =-Inf, label = r.text,hjust=1,vjust=-0.2) +
      theme_classic() + theme(text = element_text(size = 8)) + 
      theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(legend.position = 'none')
  return(p)
}

imagesc <- function(X,caxis_name='',cmap='plasma',caxis_labels=NULL,clim=c(min(X,na.rm=T),max(X,na.rm=T)),
  xlabel='',ylabel='',yticklabels=rownames(X),xticklabels=as.character(colnames(X))){
  # INPUTS:
  # X: matrix with dim names
  #
  # OUTPUTS:
  # heatmap plot of X in style of matlab imagesc

  if(is.null(caxis_labels)){ # default axis label breaks
    caxis_breaks<-labeling::extended(clim[1], clim[2], m = 5)
    caxis_labels<-as.character(labeling::extended(clim[1], clim[2], m = 5))
  } else if(is.character(caxis_labels)){ # if axis is discrete then autogenerate breaks and label with provided labels
    caxis_breaks<-labeling::extended(clim[1], clim[2], m = length(caxis_labels))
  }
  X[X>max(clim)] <- max(clim) # threshold data based on color axis
  X[X < min(clim)] <- min(clim)

  melt_mat <- melt(t(X))
  melt_mat$Var2[is.na(melt_mat$Var2)] <- 'NA'
  melt_mat$Var1 <- as.character(melt_mat$Var1)
  p<-ggplot() + geom_tile(data = melt_mat, aes(x=Var1, y=Var2, fill=value)) + 
    scale_x_discrete(limits=as.character(colnames(X)),labels=xticklabels,expand = c(0,0)) +
    scale_y_discrete(limits=rev(rownames(X)),labels=rev(yticklabels),expand = c(0,0))
  if(cmap =='plasma'){
    p <- p + scale_fill_viridis(option = 'plasma',name=caxis_name,limits=clim,breaks=caxis_breaks,labels=caxis_labels)
  } else if(cmap == 'redblue'){
    p <- p + scale_fill_gradientn(colours = c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'),
                           guide = "colorbar", limits=clim,
                           na.value = 'grey',name=caxis_name)
  } else {
    pal.idx <- which(rownames(brewer.pal.info) == cmap)  
    cols <- brewer.pal(brewer.pal.info$maxcolors[pal.idx], cmap)
    p <- p + scale_fill_gradientn(colours = cols,
                           guide = "colorbar", limits=clim,
                           na.value = 'grey',name=caxis_name,breaks=caxis_breaks,labels=caxis_labels)
  } 
  p <- p + theme_bw()
  p <- p + xlab(xlabel)+ylab(ylabel)+
      theme(text=element_text(size=8))
  return(p)

}