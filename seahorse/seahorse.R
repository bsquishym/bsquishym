# Required packages

library(TDA)
library(ggplot2)
library(ncdf4)
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

#colors for the seahorse
hmm = colorRampPalette(c("lightskyblue","firebrick1"))

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
  spot = which.min(d) #finding the minimum
  fixed = c(d[spot:length(d)],d[1:(spot-1)]) #shifting the curve to set minimums at the ends
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
  if (ln > max(d)){ #this returns a zero if we look "too high"
    return(0)
  }
  if (length(d) == 144){ # this uses fix_series if we use all longitudes
    d = fix_series(d)    # otherwise it doesn't need to
  }
  d[d < ln] = ln # setting all values less than height to the height
  # so that no odd things happen
  tmp = gridDiag(FUNvalues = d, sublevel = FALSE)[['diagram']] # Morse Filtration step
  bd = tmp[,3] - tmp[,2] # These are the overall "lifetimes"
  tmp = tmp[,2] # These are the deaths
  ft = length(tmp[tmp == ln]) # Records everything that dies on the preset spot
  ft = length(tmp[(tmp == ln) & bd > l]) # If l is used, this removes features that
  # are small than l
  return(ft) # returns feature count on a single height
}



######################################################################
#
#
#
#
#
#
######################################################################
count_gd_heights = function(d,day,lts,l,omax=F,nheight=50,lvls=NULL){
  # This counts a series of heights on a certain day.
  # d = a 3-dimensional matrix. Since our matrices that we are dealing with are
  # generally of the structure matrix[longitude,latitude,day]
  # day = the day of the year you want to look at
  # lts = the latitude to look at (from 1->n in the matrix)
  # l = the "cutoff" if we want to ignore certain small lifetimes
  # nheight = how many heights we "scan" through. default 100
  # omax is for when we do the horseZ so that every matrix that we get
  #    is on the same heights.
  # lvls is for a custome range of heights
  # and is used when creating the multiple level plots
  day_max = max(d[,,day]) # Regular max
  day_min = min(d[,,day]) # Regular min
  if (omax == T){ # omax is used when building horseZ comparison data
    # so that all matrices fall in the same format
    all_yr = paste("y",seq(1981,2021,by=1),sep="") # creates a vector that lists all years
    day_max = max(eval(parse(text = all_yr))) + 50 # finds the overall max
    day_min = min(eval(parse(text = all_yr))) - 50 # finds the overall min
  }
  if (!is.null(lvls)){ # this flag catches if we use custom heights
    day_min = lvls[1] - 50 # custom height min minus some
    day_max = lvls[2] + 50 # custom height max plus some
  }
  fts = NULL # Initializing feature vector
  hts = seq(day_min,day_max,length.out = nheight) # Initializing the heights we are looking over
  for (i in hts){ # Starting for loop
    fts = c(fts,count_gd(d[,lts,day],i,l)) # Performing count_gd on each height
    # that is in the hts vector
  }
  return(fts) # returns results
}



######################################################################
#
#
#
#
#
#
######################################################################
count_gd_day = function(d,day,plot=F,title="",l,summary=F,omax=F,nheight=50,startlat,lvls=NULL){
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
  # startlat is to keep the graphics looking right. Takes the first latitude
  # lvls is for custom range of heights
  if (omax == F){ # Flag to catch if we are using overall max or not
    day_max = max(d[,,day]) # regular max
    day_min = min(d[,,day]) # regular min
  }
  else{ # Flag that uses overall min
    all_yr = paste("y",seq(1981,2021,by=1),sep="") #creates a vector that lists all years
    day_max = max(eval(parse(text = all_yr))) + 50 # finds the overall max
    day_min = min(eval(parse(text = all_yr))) - 50 # finds the overall min
  }
  if (!is.null(lvls)){ # Flag that catches if lvls is used
    day_min = lvls[1] - 50 # custom height min minus some
    day_max = lvls[2] + 50 # custom height max plus some
  }
  n = ncol(d[,,day]) # This detects the latitudes
  mtx = NULL #Initializing the matrix that will eventually become the seahorse
  for (i in 1:n){ #this loop creates all the latitudes
    mtx = cbind(mtx,count_gd_heights(d,day,i,l,omax=omax,nheight=nheight,lvls=lvls))
    # This is the part that loops through and creates all of the features by repeatedly
    # calling count_gd_heights for each latitude
    if (summary == T){ # Flag that catches if we should overlay any information
      sum_stats = apply(d,2,summary) # summary stats for each latitude
      mns = sum_stats[4,] # the means of each latitude
      sds = apply(d,2,sd) # the standard deviations of each latitude
      mdns = sum_stats[3,] # the medians of each latitude
      q1 = sum_stats[2,]# first quantile
      q3 = sum_stats[5,] # third quantile
    }
  }
  if (plot == T){ #plots the information
    mtxp = mtx # creates a copy of the matrix
    mtxp[mtxp >= 5] = 5 #fixes color range
    hts = seq(day_min,day_max,length.out = 15) # the heights to plot
    tixy = round(hts,0) # labels that are rounded values so it's not ugly
    #
    # A little explanation here
    # since the image.plot function from fields forces everything to be plotted on a [0,1]
    # interval, we have to scale everything to that.
    # this small process normalizes all the labels and such so the axis labels work
    # and look alright
    tixy_at = (hts - min(hts))/(max(hts) - min(hts)) # normalizing y-axis value locations
    tixx = seq(startlat,startlat+(n-1)*2.5,by=2.5) # labels for x-axis
    tixx_at = (tixx - min(tixx))/((startlat+(n-1)*2.5) - min(tixx)) # normalizing x-axis value locations
    colrange = range(mtx)[2] - range(mtx)[1] # choosing the range of colors
    #
    # This is the part that plots the seahorse
    # We use the transpose so that we get the wanted orientation. Otherwise it's more vertical and
    # looks more like a sprinkler or waterspout
    # This part could be modified to display the information using whatever graphical package
    # that you wanted
    # Important parts of this are:
    #  The breaks and axis.args are used to make the legend more readable
    #  axes = FALSE gets rid of original axes so we can overlay ours
    #  cex.main = 2 makes the overall title larger because its difficult to see sometimes
    #  the two axis lines at the end are the ones that overlay our custom axis
    image.plot(t(mtxp),col=c("white",hmm(4),"firebrick4"),
               breaks = c(-1,0,1,2,3,4,5),
               axes=FALSE,
               ylab = "Geopotential Height", xlab = "Latitude",
               main=title,
               xlim = c(0,1),ylim=c(0,1),
               axis.args = list(at = seq(-0.5,4.5,by=1), labels = c(0,1,2,3,4,">=5")),
               cex.main = 2)
    axis(side = 1, at = tixx_at, labels = tixx) # custom x axis
    axis(side = 2, at = tixy_at, labels = tixy) # custom y axis
    if (summary == T){ # Flag catches if we want the summary overlaid
      # I'll be honest I haven't used this much. Hopefully nothing breaks in it
      sdd_f = cbind(mns - 2*sds, mns - sds, mns + sds, mns + 2*sds) # band of 2 sds
      mns = (mns - day_min+l)/(day_max - day_min+2*l) # normalizes means
      sdd_f = (sdd_f - day_min+l)/(day_max - day_min+2*l) #normalizes the above band
      mdns = (mdns - day_min)/(day_max - day_min) #normalizes medians
      q1 = (q1 - day_min)/(day_max - day_min) #normalizes first quantile
      q3 = (q3 - day_min)/(day_max - day_min) #normalizes third quantile
      tst = cbind(mns,tixx_at) # binds to above custom axis
      sdd = cbind(sdd_f,tixx_at) #binds to above custom axis
      tst2 = cbind(mdns,tixx_at) #binds to above custom axis
      qq1 = cbind(q1,tixx_at) #binds to above custom axis
      qq3 = cbind(q3,tixx_at) #binds to above custom axis
      lines(tst[,2],tst[,1],lwd=2,lty=2,col="black") # all of these line comments
      lines(sdd[,5],sdd[,1],lwd=2,lty=3,col="black") # add lines to the plot
      lines(sdd[,5],sdd[,2],lwd=2,lty=3,col="black")
      lines(sdd[,5],sdd[,3],lwd=2,lty=3,col="black")
      lines(sdd[,5],sdd[,4],lwd=2,lty=3,col="black")
      lines(tst2[,2],tst2[,1],lwd=2,lty=6,col="black")
      lines(qq1[,2],qq1[,1],lwd=2,lty=3,col="darkred")
      lines(qq3[,2],qq3[,1],lwd=2,lty=3,col="darkred")
    }
  }
  return(mtx) # returns the matrix that the seahorse is plotted from
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
  
  # The y3000 part is something I broke that I don't recall how to fix so
  # this is like a weird bandaid
  y3000 = "don't worry about this"
  if (years_back != 0){ # This flag catches to correctly create a sequence of years
    yrs_scan = seq(year_comp,year_comp - years_back,by = -1) # sequence of years
  }
  else{ # this is part of the bandaid above
    yrs_scan = 3000 # part of the bandaid
  }
  if (daycomp - days_back <= 0){ # This catches if you look at a range of days that
    # goes into a previous year
    day_scan = seq(daycomp-1,1,by = -1) # sequence of days to scan
    day_scan_other = seq(365,365 + daycomp - days_back,by = -1) # the second part that is needed
    # since it overlaps into previous year
  }
  else{
    day_scan = seq(daycomp,daycomp-days_back,by = -1) # range of days if everything is fine
  }
  backyrs = NULL # initializing the years that will be created
  
  # (***)
  # This next loop is slightly strange maybe
  # It goes through the years in yrs_scan and gathers all the days into a single
  # 3-Dimensional matrix. It uses abind to "stack" them up.
  # Eventually creating a single matrix that has every piece of data needed
  # This is done because of how I set up the previous functions and can likely
  # be changed to be more computationally efficient
  for (j in yrs_scan){ # loop initialize
    tmp_hold = eval(parse(text = paste("y",j,sep=""))) # temporarily holds the year
    if (j == year_comp & j != y3000){ # part of the bandaid above
      backyrs = abind(backyrs,tmp_hold[,,day_scan],along=3) # stacking days
    }
    else{
      backyrs = abind(backyrs,tmp_hold[,,day_scan],along=3) # stacking days
      if (daycomp - days_back <= 0){
        backyrs = abind(backyrs,tmp_hold[,,day_scan_other],along=3) # stacking days if overlapping
      }
    }
  }
  # This if statement catches if we overlap years
  if (daycomp - days_back <= 0){ # seeing if we overlap years
    last_yr = min(yrs_scan) - 1 # looking at an additional year since we overlap
    tmp_hold = eval(parse(text = paste("y",last_yr,sep="")))
    backyrs = abind(backyrs,tmp_hold[,,day_scan_other],along=3)
  }
  mtx = NULL # initializing the empty matrix to put everything
  n = dim(backyrs)[3] # This is the number of TOTAL DAYS that are being looked at
  # or the number of "stacked" days from (***) above
  for (i in 1:n){ # looping and calling count_gd_days for each day in the stack
    mtx = abind(mtx,t(count_gd_day(backyrs,i,plot=F,l=l,omax=T,startlat=0)),along=3)
  }
  
  # So we know that apply(something,1,blah) computes blah on each row
  # and apply(something,2,blah) computes blah on each column
  # well apply(something,c(1,2),blah) computes blah no each cell through the stack
  # For example, if we have something that is dimension 3x3x9, the apply will start at
  # (1,1,1) location and go through to (1,1,9) and compute the blah for that particular "cell"
  # in our case, mean and sd
  means = apply(mtx,c(1,2),mean) #mean of each cell
  sds = apply(mtx,c(1,2),sd) #sd of each cell
  fin = list("Mean" = means, "SD" = sds,"Raw" = mtx) #putting this all in a list
  return(fin) # returning the list
}



######################################################################
#
#
#
#
#
#
######################################################################
horsieZ = function(d,daycomp,mtx,l,plotting=F,title="",startlat,lvls=NULL){
  # This compares a particular day to the base created above
  # d = dataset of the year
  # daycomp = day to compare
  # mtx = the "base" that you created above
  # l -> set at 0 *broken*
  # plotting = if you want to plot. returns only matrix otherwise
  # title is... the title of the plot.
  # startlat is to plot it correctly if we use different latitudes
  # lvls is for custom heights (similar to other functions)
  hmmNEG = colorRampPalette(c("navyblue","lightskyblue")) #colors for negative z
  hmmPOS = colorRampPalette(c("khaki1","firebrick1")) #colors for positive z
  all_yr = paste("y",seq(1981,2021,by=1),sep="") # creating vector of all years
  day_max = max(eval(parse(text = all_yr))) + 50 # finding overall max
  day_min = min(eval(parse(text = all_yr))) - 50 # finding overall mean
  if (!is.null(lvls)){ # Flag that catches for custom heights
    day_min = lvls[1] - 50 #minimum custom minus some
    day_max = lvls[2] + 50 #maximum custom plus some
  }
  means = mtx$Mean # extracting means from base data list
  sds = mtx$SD # extracting standard deviations from base data list
  n = ncol(d[,,daycomp]) # number of latitudes
  
  # The following two statements create the "comparison" seahorse
  if (is.null(lvls)){ # checks if custom heights are used
    mtx_comp = t(count_gd_day(d,daycomp,plot=F,l=l,omax=T))
  }
  else{ # this is used if custom heights are needed
    mtx_comp = t(count_gd_day(d,daycomp,plot=F,l=l,omax=F,lvls=lvls))
  }
  mtx_z = (mtx_comp - means)/sds # calculates Z-score for each cell
  
  # These computations are fixing several cases where the Z-score is infinite
  # or divided by zero or something of that nature
  # it also takes very large values and places them into a single value
  # We do this to make the coloring on the plot more  uniform
  # omitting this part gives you greater control of seeing extremely unlikely events
  # with horrifically large z-scores, but makes the coloring very wonky
  # due to the greater range
  mtx_z[is.infinite(mtx_z)] = 4.1
  mtx_z[is.nan(mtx_z)] = NA
  mtx_z[mtx_z > 4] = 4.1
  mtx_z[mtx_z < -4] = -4.1
  if (plotting == T){ # Flag to check if we are plotting
    tixy = seq(round(day_min+l,-2),round(day_max+l,-2),by=100) # y-axis label values
    tixy_at = (tixy - min(tixy))/(max(tixy) - min(tixy)) # y-axis location
    tixx = seq(startlat,startlat+(n-1)*2.5,by=2.5) # x-axis label values
    tixx_at = (tixx - min(tixx))/((startlat+(n-1)*2.5) - min(tixx)) # x-axis locations
    colrange = range(mtx)[2] - range(mtx)[1] # range of colors
    
    #
    #
    # Similar to count_gd_day above
    # Except the matrix is already transposed so we don't need it
    # This part could be modified to display the information using whatever graphical package
    # that you wanted
    # Important parts of this are:
    #  The breaks and axis.args are used to make the legend more readable
    #  axes = FALSE gets rid of original axes so we can overlay ours
    #  cex.main = 2 makes the overall title larger because its difficult to see sometimes
    #  the two axis lines at the end are the ones that overlay our custom axis
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
  fin = list("Mean" = means, "SD" = sds,"Z" = mtx_z,"Raw" = mtx) # creating a list
  return(mtx_z) # returns a list if needed
}


######################################################################
#
#
#
#
#
#
######################################################################
build_entire_comp = function(){
  #
  # This function builds a comparison database with whatever years you have loaded
  # years must be named "yXXXX"
  # For example, 2020 should be named "y2020"
  # This does NOT take into account leap years
  # This does NOT take into account leap years
  # This does NOT take into account leap years
  # This does NOT take into account leap years
  #
  # We do not account for leap years as I am slightly lazy apparently
  #
  # It returns a very large list. Let data = build_entire_comp()
  #   Then data[[1]] is January
  #   data[[2]] is February
  #   ...
  #   data[[12]] is December
  #   Each of these is the same as build_comp from above
  #   So data[[1]]$Raw, data[[1]]$Mean, and data[[1]]$SD exist
  #
  # This can take a long time to run, hence the print statements to see where you are
  # For example, when doing the entire northern hemisphere with nheight=50
  # It takes about 1.25 hours
  # This is due to both R being slow, and probably me being inefficient at coding this at the
  # start of research and not going back and fixing it
  # This can very likely be improved drastically
  start = Sys.time() # This is for tracking the amount of time it takes
  print("Building January...")
  jan30_0 = build_comp(daycomp=31,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building February...")
  feb30_0 = build_comp(daycomp=59,year_comp=2010,days_back=27,years_back=29,l=0)
  print("Building March...")
  mar30_0 = build_comp(daycomp=90,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building April...")
  apr30_0 = build_comp(daycomp=120,year_comp=2010,days_back=29,years_back=29,l=0)
  print("Building May...")
  may30_0 = build_comp(daycomp=151,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building June...")
  jun30_0 = build_comp(daycomp=181,year_comp=2010,days_back=29,years_back=29,l=0)
  print("Building July...")
  jul30_0 = build_comp(daycomp=212,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building August...")
  aug30_0 = build_comp(daycomp=243,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building September...")
  sep30_0 = build_comp(daycomp=273,year_comp=2010,days_back=29,years_back=29,l=0)
  print("Building October...")
  oct30_0 = build_comp(daycomp=304,year_comp=2010,days_back=30,years_back=29,l=0)
  print("Building November...")
  nov30_0 = build_comp(daycomp=334,year_comp=2010,days_back=29,years_back=29,l=0)
  print("Building December...")
  dec30_0 = build_comp(daycomp=365,year_comp=2010,days_back=30,years_back=29,l=0)
  fin = list(jan30_0,feb30_0,mar30_0,apr30_0,
             may30_0,jun30_0,jul30_0,aug30_0,
             sep30_0,oct30_0,nov30_0,dec30_0) # creating a list
  end = Sys.time() # for tracking time
  runtime = end - start # computing total runtime
  print(runtime) # outputting the runtime
  return(fin) # returning a giant list
}

