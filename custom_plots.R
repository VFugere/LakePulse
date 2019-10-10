mybxp <- function(x,y,colvec,ylab){
  boxplot(y~x, lty=1, staplewex=0, whisklwd=1,boxwex=0.6, outline=F,boxlwd=1, medlwd=1, col=colvec, ylab=ylab)
  points(x=jitter(as.numeric(x)),y=y, pch=16,cex=0.5, col=alpha(1,0.3))
}

mybubble <- function(x,y,z){
  plot(y~x, bty='n', xlab='Lake area',ylab='human impact index',type='n')
  cexvec <- scales::rescale(z, to=c(0.5,3.5))
  points(y~x, cex=cexvec, pch=16, col = scales::alpha('navy blue',0.5))
}
