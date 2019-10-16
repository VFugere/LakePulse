mybxp <- function(x,y,colvec,ylab){
  boxplot(y~x, lty=1, staplewex=0, whisklwd=1,boxwex=0.6, outline=F,boxlwd=1, medlwd=1, col=colvec, ylab=ylab)
  points(x=jitter(as.numeric(x)),y=y, pch=16,cex=0.5, col=alpha(1,0.3))
}

mybubble <- function(x,y,z,name='',ez=F){
  plot(y~x, bty='l', ylab='human impact index',ylim=c(0,1),xlab = 'lake area',type='n',log='x',xaxt='n',yaxt='n')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0.1,1,10,100),lab=c('0.1','1','10','100'))
  cexvec <- scales::rescale(z, to=c(0.5,3.5))
  if(name != ''){
  text(x=10^1.5,y=0.9,pos=3,label=make.italic(name))
  points(x=c(10^1,10^1.25,10^1.5,10^1.75,10^2),y=rep(0.85,5),pch=1,cex=seq(0.5,3.5,length.out = 5))
  }
  if(length(ez) == 1){
  points(y~x, cex=cexvec, pch=16, col = scales::alpha('navy blue',0.5))
  }else{
  points(y~x, cex=cexvec, pch=16, col = scales::alpha(cols2[as.numeric(ez)],0.5))
  legend('right',bty='n',legend=levels(ez),pch=15,col=cols2)
  }
}

mybubble2 <- function(x,y,z,name='',ez=F){
  plot(y~x, bty='l', ylab='human impact index',ylim=c(0,1),type='n',log='x',xaxt='n',yaxt='n')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  cexvec <- scales::rescale(z, to=c(0.5,3.5))
  if(name != ''){
    legend('topright',bty='n',legend=make.italic(name))
    #points(x=c(10^1,10^1.25,10^1.5,10^1.75,10^2),y=rep(0.85,5),pch=1,cex=seq(0.5,3.5,length.out = 5))
  }
  if(length(ez) == 1){
    points(y~x, cex=cexvec, pch=16, col = scales::alpha('navy blue',0.5))
  }else{
    points(y~x, cex=cexvec, pch=16, col = scales::alpha(cols2[as.numeric(ez)],0.5))
    legend('right',bty='n',legend=levels(ez),pch=15,col=cols2)
  }
}

mybubble3 <- function(x,y,z,name='',ez=F){
  y <- y * 100
  plot(y~x, bty='l', xlab='human impact index',ylab = name,type='n',xaxt='n',yaxt='n',xlim=c(0,1))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  cexvec <- scales::rescale(z, to=c(0.5,3.5))
  points(y~x, cex=cexvec, pch=16, col = scales::alpha(cols2[as.numeric(ez)],0.5))
  legend('right',bty='n',legend=levels(ez),pch=15,col=cols2)
  text(x=0.875,y=max(y)*0.95,pos=3,label=make.italic('lake depth'))
  points(x=seq(from=0.75,to=1,length.out = 5),y=rep(max(y)*0.9,5),pch=1,cex=seq(0.5,3.5,length.out = 5))
}

mapplot <- function(x,y,z,name){
  
  xrange <- range(x)+c(-2,2)
  yrange <- range(y)+c(-1,1)
  plot(map, xlim = xrange, ylim = yrange,col='light gray',border=0,asp=1.2,axes=F,cex.lab=0.5)
  
  colfunc <- colorRampPalette(RColorBrewer::brewer.pal(11,'RdYlBu'))
  mapcols <- colfunc(100)[100:1]
  zvec <- scales::rescale(z, to=c(0,100))
  zvec <- round(zvec,0)
  
  points(x=x,y=y,pch=1,col=1,cex=1.2)
  points(x=x,y=y,pch=16,col=alpha(mapcols[zvec],0.8),cex=1.1)
  
  legend('topright',legend=make.italic(name),bty='n',text.col=1)
  
  leg.x.rg <- range(xrange)[1]-(range(xrange)[1]-range(xrange)[2])*0.75
  leg.x.rg <- seq(leg.x.rg,range(xrange)[2],length.out = 100)
  ypos <- range(yrange)[1]
  points(rep(ypos[1],100)~leg.x.rg, col=mapcols, pch=16)
  #text(x=leg.x.rg[50],y=ypos+0.2,cex=1,label=make.italic(name),pos=3)
  text(x=leg.x.rg[c(1,100)],y=ypos-0.2,cex=1,label=c('low','high'),pos=1)

}


