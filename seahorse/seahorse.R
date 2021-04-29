# Required packages

library(TDA)
library(ncdf4)
library(tidyverse)
library(plotly)
library(abind)
library(fields)

#
#
# included is a file called bymonth30.R
# This is a 30-year horsieZ base with l=50
# to load it, use load(file="bymonth30.RData")
# it will import it as variable name bymonth30
#
# you can see specific months via bymonth30[[i]] where i = 1:12
# where the months are numbered as you would expect
# additionally, bymonth30[[i]]$Mean, bymonth30[[i]]$SD, and
#   bymonth30[[i]]$Raw exist for all months
# you can create your own bases but this can save computation time
#   as creating this base took a bit of time, which can increase if
#   you increase nheight to make a finer mesh
#
#
#
# There is currently erroneaously display of geopotential heights
# when using l is large. The data is correct but it is difficult to display the locations
# correctly as certain features will eventually having a larger life than l, but will not be counted
# until l is reached. This can make it appear that features aren't born until they reach "maturity" I guess?
# anyway I'm unsure how to fix it at the moment but I will come back to it.
# Thaaaaaaanks!

######################################################################
#
#
#
#
#
#
######################################################################
fix_series = function(d){
  # this small function fixes the series so that the minimum is on the
  # borders. This is to make it so that the features are not accidentally
  # counted twice.
  
  # it takes an input of a continuous one dimensional series/curve
  spot = which.min(d)     
  fixed = c(d[spot:length(d)],d[1:(spot-1)])
  return(fixed)
}




######################################################################
#
#
#
#
#
#
######################################################################
count_gd = function(d,ln,l){
  # This counts the number of features when we set the height
  # d = a series (one latitude on one day, or one dimensional curve)
  # ln = the height that we are checking
  # l = the "cutoff" if we want to ignore certain small lifetimes
  if (ln > max(d)){
    return(0)
  }
  d = fix_series(d)
  d[d < ln] = ln
  tmp = gridDiag(FUNvalues = d, sublevel = FALSE)[['diagram']]
  bd = tmp[,3] - tmp[,2]
  tmp = tmp[,2]
  ft = length(tmp[tmp == ln])
  ft = length(tmp[(tmp == ln) & bd > l])
  return(ft)
}



######################################################################
#
#
#
#
#
#
######################################################################
count_gd_heights = function(d,day,lts,l,omax=F,nheight=100){
  # This counts a series of heights on a certain day.
  # d = a 3-dimensional matrix. Since our matrices that we are dealing with are
  # generally of the structure matrix[longitude,latitude,day]
  # day = the day of the year you want to look at
  # lts = the latitude to look at (from 1->n in the matrix)
  # l = the "cutoff" if we want to ignore certain small lifetimes
  # nheight = how many heights we "scan" through. default 100
  # omax is for when we do the horseZ so that every matrix that we get
  #    is on the same heights.
  day_max = max(d[,,day])
  day_min = min(d[,,day])
  if (omax == T){
    day_max = 6018
    day_min = 4594
  }
  fts = NULL
  hts = seq(day_min,day_max,length.out = nheight)
  # hts is which heights to check
  for (i in hts){
    fts = c(fts,count_gd(d[,lts,day],i,l))
  }
  return(fts)
}



######################################################################
#
#
#
#
#
#
######################################################################
count_gd_day = function(d,day,plot=F,title="",l,summary=F,omax=F,nheight=100){
  # This is the main function that uses the others
  # d = 3-dimensional matrix of data
  # day = day that will be looked at
  # plot = If the seahorse plot will be added or not
  # title = Title on the seahorse plot
  # l = the "cutoff" if we want to ignore certain small lifetimes
  # summary = If means/sds/medians will be plotted
  # nheight = how many heights we "scan" through. default 100
  # omax is for when we do the horseZ so that every matrix that we get
  #    is on the same heights.
  hmm = colorRampPalette(c("lightskyblue","firebrick1")) #colors
  
  if (omax == F){
    day_max = max(d[,,day])
    day_min = min(d[,,day])
  }
  else{
    day_max = 6018
    day_min = 4594
  }
  mns = NULL
  mdns = NULL #various stats to overlay
  q1 = NULL
  q3 = NULL
  sds = NULL
  n = ncol(d[,,day])
  mtx = NULL
  for (i in 1:n){ #this loop creates all the latitudes
    mtx = cbind(mtx,count_gd_heights(d,day,i,l,omax=omax,nheight=nheight))
    if (summary == T){
      mns = c(mns,mean(d[,i,day]))
      sds = c(sds,sd(d[,i,day]))
      mdns = c(mdns,median(d[,i,day]))
      q1 = c(q1,as.numeric(quantile(d[,i,day])[2]))
      q3 = c(q3,as.numeric(quantile(d[,i,day])[4]))
    }
  }
  if (plot == T){ #plots the information
    mtxp = mtx
    mtxp[mtxp >= 5] = 5 #fixes color range
    hts = seq(day_min,day_max,length.out = 15)
    tixy = round(hts,0)
    print(tixy)
    tixy_at = (hts - min(hts))/(max(hts) - min(hts))
    print(tixy_at)
    tixx = seq(42.5,80,by=2.5)
    tixx_at = (tixx - min(tixx))/(80 - min(tixx))
    colrange = range(mtx)[2] - range(mtx)[1]
    image.plot(t(mtxp),col=c("white",hmm(4),"firebrick4"),
               breaks = c(-1,0,1,2,3,4,5),
               axes=FALSE,
               ylab = "Geopotential Height", xlab = "Latitude",
               main=title,
               xlim = c(0,1),ylim=c(0,1),
               axis.args = list(at = seq(-0.5,4.5,by=1), labels = c(0,1,2,3,4,">=5")),
               cex.main = 2)
    axis(side = 1, at = tixx_at, labels = tixx)
    axis(side = 2, at = tixy_at, labels = tixy)
    if (summary == T){
      sdd_f = cbind(mns - 2*sds, mns - sds, mns + sds, mns + 2*sds)
      mns = (mns - day_min+l)/(day_max - day_min+2*l)
      sdd_f = (sdd_f - day_min+l)/(day_max - day_min+2*l)
      mdns = (mdns - day_min)/(day_max - day_min)
      q1 = (q1 - day_min)/(day_max - day_min)
      q3 = (q3 - day_min)/(day_max - day_min)
      tst = cbind(mns,tixx_at)
      sdd = cbind(sdd_f,tixx_at)
      tst2 = cbind(mdns,tixx_at)
      qq1 = cbind(q1,tixx_at)
      qq3 = cbind(q3,tixx_at)
      lines(tst[,2],tst[,1],lwd=2,lty=2,col="black")
      lines(sdd[,5],sdd[,1],lwd=2,lty=3,col="black")
      lines(sdd[,5],sdd[,2],lwd=2,lty=3,col="black")
      lines(sdd[,5],sdd[,3],lwd=2,lty=3,col="black")
      lines(sdd[,5],sdd[,4],lwd=2,lty=3,col="black")
      lines(tst2[,2],tst2[,1],lwd=2,lty=6,col="black")
      lines(qq1[,2],qq1[,1],lwd=2,lty=3,col="darkred")
      lines(qq3[,2],qq3[,1],lwd=2,lty=3,col="darkred")
    }
  }
  return(mtx)
}



######################################################################
#
#
#
#
#
#
######################################################################
build_comp = function(daycomp,year_comp,days_back,years_back,l){
  # This creates the "base" that we compare against
  # It assumes that years are given variable name yXXXX where XXXX is the year.
  #   For example, 2010 would be y2010
  #
  #
  # daycomp = The first day to "start at"
  # year_comp = The first year to "start at"
  # days_back = How far to go back.
  #   For example, daycomp=31, days_back=30 will do days 1-31 (January)
  # years_back = Same as day but for year
  # l = should be set to zero.
  # This also accounts for going past years. Such as daycomp=2 days_back=5
  #   Won't break so yay!
  #
  # Outputs a list
  # $Mean -> Mean of each row/column entry
  #   For example, if the matrix is 100x100 we will have 10,000 means
  # $SD -> SD for each row/column entry
  # $Raw -> the raw data that computes the mean and sd. Could be useful
  #
  
  y3000 = "don't worry about this"
  if (years_back != 0){
    yrs_scan = seq(year_comp,year_comp - years_back,by = -1)
  }
  else{
    yrs_scan = 3000
  }
  if (daycomp - days_back <= 0){
    day_scan = seq(daycomp-1,1,by = -1)
    day_scan_other = seq(365,365 + daycomp - days_back,by = -1)
  }
  else{
    day_scan = seq(daycomp,daycomp-days_back,by = -1)
  }
  backyrs = NULL
  for (j in yrs_scan){
    print(j)
    tmp_hold = eval(parse(text = paste("y",j,sep="")))
    if (j == year_comp & j != y3000){
      backyrs = abind(backyrs,tmp_hold[,,day_scan],along=3)
    }
    else{
      backyrs = abind(backyrs,tmp_hold[,,day_scan],along=3)
      if (daycomp - days_back <= 0){
        backyrs = abind(backyrs,tmp_hold[,,day_scan_other],along=3)
      }
    }
  }
  if (daycomp - days_back <= 0){
    last_yr = min(yrs_scan) - 1
    tmp_hold = eval(parse(text = paste("y",last_yr,sep="")))
    backyrs = abind(backyrs,tmp_hold[,,day_scan_other],along=3)
  }
  mtx = NULL
  n = dim(backyrs)[3]
  for (i in 1:n){
    print(i) #tracking how long you have to wait
    mtx = abind(mtx,t(count_gd_day(backyrs,i,plot=F,l=l,omax=T)),along=3)
  }
  means = apply(mtx,c(1,2),mean)
  sds = apply(mtx,c(1,2),sd)
  fin = list("Mean" = means, "SD" = sds,"Raw" = mtx)
  return(fin)
}



######################################################################
#
#
#
#
#
#
######################################################################
horsieZ = function(d,daycomp,mtx,l,plotting=F,title=""){
  # This compares a particular day to the base created above
  # d = dataset of the year
  # daycomp = day to compare
  # mtx = the "base" that you created above
  # l -> set at 0 *broken*
  # plotting = if you want to plot. returns only matrix otherwise
  # title is... the title of the plot.
  hmmNEG = colorRampPalette(c("navyblue","lightskyblue")) #colors for negative z
  hmmPOS = colorRampPalette(c("khaki1","firebrick1")) #colors for positive z
  day_min = 4594
  day_max = 6018
  means = mtx$Mean
  sds = mtx$SD
  mtx_comp = t(count_gd_day(d,daycomp,plot=F,l=l,omax=T))
  mtx_z = (mtx_comp - means)/sds
  mtx_z[is.infinite(mtx_z)] = 4.1
  mtx_z[is.nan(mtx_z)] = NA
  mtx_z[mtx_z > 4] = 4.1
  mtx_z[mtx_z < -4] = -4.1
  if (plotting == T){
    tixy = seq(round(day_min+l,-2),round(day_max+l,-2),by=100)
    tixy_at = (tixy - min(tixy))/(max(tixy) - min(tixy))
    tixx = seq(42.5,80,by=2.5)
    tixx_at = (tixx - min(tixx))/(80 - min(tixx))
    colrange = range(mtx)[2] - range(mtx)[1]
    image.plot(mtx_z,col=c(hmmNEG(4),"lightsteelblue1","oldlace",hmmPOS(4)),
               breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
               axes=FALSE,
               ylab = "Geopotential Height", xlab = "Latitude",
               main=title,
               xlim = c(0,1),ylim=c(0,1),
               axis.args = list(at = seq(-4.5,4.5,by=1),
                                labels = c("< -4","(-4 , -3)","(-3, -2)","(-2 , -1)",
                                           "(-1 , 0)","(0 , 1)",
                                           "(1 , 2)","(2 , 3)","(3 , 4)","> 4")),
               cex.main=2)
    axis(side = 1, at = tixx_at, labels = tixx)
    axis(side = 2, at = tixy_at, labels = tixy)
  }
  # image.plot(mtx[,,1])
  fin = list("Mean" = means, "SD" = sds,"Z" = mtx_z,"Raw" = mtx)
  return(mtx_z)
}


######################################################################
#
#
#
#
#
#
######################################################################
