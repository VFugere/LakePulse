mybxp <- function(x,y,colvec,ylab){
  boxplot(y~x, lty=1, staplewex=0, whisklwd=1,boxwex=0.6, outline=F,boxlwd=1, medlwd=1, col=colvec, ylab=ylab)
  points(x=jitter(as.numeric(x)),y=y, pch=16,cex=0.5, col=alpha(1,0.3))
}