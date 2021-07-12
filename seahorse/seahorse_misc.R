# Requires having all the information in seahorse.R active
# Random Comment

######################################################################
#
# qk_series will output a series over a particular range of days
# this allows us to look at a time series over any range of time
# for any computable summary statistic
#
#
######################################################################

qk_series = function(d,l=0,func){
  # d = the year we want to look at, as well as the range of days
  #     for example, y1990[,,1:30] will look at days 1-30
  # l = DO NOT USE
  # func = the function we want to use to evaluate
  #   for example, max will take the max features of the day
  #   mean will take the mean
  #   median will take the median
  #   etc...
  srs = NULL # Initializing series
  # This small if statement checks if we get a matrix without a 3rd dimension
  # AKA it's 30x30 instead of 30x30x5 or something
  # If it's not, it forces it to be one by repeating itself
  if (is.na(dim(d)[3])){
    d = abind(d,d,along=3)
    srs = count_gd_day(d,1,plot=F,l=0)
  }
  else{
    n = dim(d)[3] #the days of the 3-dimensional dataset
    for (i in 1:n){
      tmp = count_gd_day(d,i,plot=F,l=0,startlat=0) # calling count_gd_day
      # we can use startlat = 0
      srs = c(srs,func(tmp)) # creating vector that uses function on matrix from above
    }
  }
  return(srs)
  # This function will be used in basically all of the series
  # that we create
}



######################################################################
#
# count_gd_season performs qk_series over all seasons
# We define seasons as
#   winter -> dec-feb
#   spring -> mar-may
#   sum -> jun-aug
#   fall -> sep-nov
#
# If we choose the year 2010, it takes december of 2009, and
#   the rest from 2010
######################################################################
count_gd_season = function(d,l=0,func,norm=F){
  # d = the year we want to look at.
  # l = DO NOT USE
  # func = the function we want to use to evaluate
  #   for example, max will take the max features of the day
  #   mean will take the mean (over the non-zero feature counts)
  #   median will take the median over the non-zero feature counts
  # norm = If we want to normalize the information or not. It is mostly
  #   useless for this function
  ref_yr = str_remove(deparse(substitute(d)),"y") # reference year
  wint_yr_ind = as.numeric(ref_yr) - 1 # winter year
  wint_yr = eval(parse(text = paste("y",wint_yr_ind,sep=""))) # getting winter year
  if (dim(wint_yr)[3] == 366){ # This checks for leap year and slices accordingly
    wint_yr = wint_yr[,,336:366]
    d_wint = d[,,1:59]
    d_spring = d[,,60:151]
    d_summer = d[,,152:244]
    d_fall = d[,,245:334]
  }
  else{ # This runs if there is no leap year
    wint_yr = wint_yr[,,335:365]
  }
  if (norm == T){ # This is normalizing over the year
    dmax = max(d) # getting max
    dmin = min(d) # getting min
    d = (d - dmin)/(dmax - dmin) # normalizing
    wint_yr = (wint_yr - dmin)/(dmax - dmin) # normalizing the winter year with current
    # year's min/max
  }
  if (dim(d)[3] == 366){ # More leap year fun
    d_wint = d[,,1:60]
    d_spring = d[,,61:152]
    d_summer = d[,,153:245]
    d_fall = d[,,246:335]
  }
  else{
    d_wint = d[,,1:59]
    d_spring = d[,,60:151]
    d_summer = d[,,152:244]
    d_fall = d[,,245:334]
  }
  d_winter = abind(wint_yr,d_wint,along=3) # creating one winter
  print("Calculating winter...")
  s_winter = qk_series(d_winter,l=0,func=func) # Creating series for the seasons
  print("Calculating spring...")
  s_spring = qk_series(d_spring,l=0,func=func)
  print("Calculating summer...")
  s_summer = qk_series(d_summer,l=0,func=func)
  print("Calculating fall...")
  s_fall = qk_series(d_fall,l=0,func=func)
  fin = list(s_winter,s_spring,s_summer,s_fall) # Creating list to return
  return(fin) #returns list of series run through the functions
}




######################################################################
#
# This looks at the lifetimes of features.
# Lifetime is thought of as the time it takes for a feature to die
#   after it is born.
# This function then sums those features per latitude and returns
#   a vector of results
#
######################################################################
life_time = function(d,day,l=0,norm=F){
  # norm = Normalizes over the day and 30-yr climatology
  #        Makes comparisons easier
  # d = 3-dimensional matrix
  # day = day of the year
  # set l = 0
  n = dim(d)[2] #the latitudes to use
  life_avg = NULL # initializing vector
  for (i in 1:n){ # looping over number of latitudes
    m = d[,i,day] # Creating temporary matrix
    if (norm == T){ # normalizing temporary matrix if needed
      omax = max(d)
      omin = min(d)
      m = (m - omin)/(omax - omin)
    }
    tmp = gridDiag(FUNvalues=fix_series(m),sublevel=FALSE)[['diagram']] # Morse filtration
    lif = tmp[,3] - tmp[,2] # lifetimes
    lif = lif[lif > l] # lifetime greater than cutoff
    life_avg = c(life_avg,sum(lif)) # currently returns sum of life
    # despite being named life_avg
    # I realize this is confusing I changed things a bit
  }
  return(life_avg) # returning the vector
}



######################################################################
#
# life_series_mnth does life_time over a series of days and performs
#   whatever function to return a vector of results
# ex if you choose mean you will get the mean of the lifetime sums
#   over what days you calculated over
#
######################################################################
life_series_mnth = function(d,days,l,func,norm=F){
  # d = 3-dimensional matrix
  # l = lifetime to ignore
  # func = function we pass. mean, max, etc...
  lives = NULL # initializing lives
  if (is.na(dim(d)[3])){ #checks if 3d matrix or not and fixes it by stacking if it isn't
    d = abind(d,d,along=3) #stacking matrix
    lives = life_time(d,1,l) #performs life_time over the single day
  }
  else{
    n = dim(d)[3] #number of days in the year
    for (i in days){
      lives = cbind(lives,life_time(d,i,l,norm=norm)) #performs life_time over series of days
    }
  }
  if (is.null(dim(lives))){ # This checks if a matrix or not. vectors return null
    fin = lives
  }
  else{ # This applies the function over the matrix to turn it into a vector
    fin = apply(lives,1,func)
  }
  return(fin)
  # this returns the function applied over the latitudes
  # lives will originally be a day x latitude matrix
  # mean will take the mean over each latitude over however many days
  # so the first entry of the output is the average lifetime
  # for the first latitude over however many days
}




######################################################################
#
# life_season does the same thing for lifetimes as it the other seasonal
#   function does for features
#
#
#
######################################################################
life_season = function(d,l,func,norm=F){
  # d = 3-dimensional matrix
  # l = the "cutoff" if we want to ignore certain small lifetimes
  # func = function we pass. mean, max, etc...
  # NOTE: if using y2010, for example, y2009 must also exist
  # since it grabs winter from that year
  #
  # norm sets if we normalize the year first
  # if norm = T, then l should be small
  # around 0.01 - 0.15
  seasons = c("Winter","Spring","Summer","Fall") # Good ol' seasons
  ref_yr = str_remove(deparse(substitute(d)),"y") # reference year
  wint_yr_ind = as.numeric(ref_yr) - 1 # winter year
  wint_yr = eval(parse(text = paste("y",wint_yr_ind,sep=""))) # importing winter year
  if (dim(wint_yr)[3] == 366){ # checking for leap year and slicing accordingly
    wint_days_prev = 336:366
    wint_days = 1:59
    spring_days = 60:151
    summer_days = 152:244
    fall_days = 245:334
  }
  else{
    wint_days_prev = 335:365
  }
  if (dim(d)[3] == 366){ # leap year fun
    wint_days = 1:60
    spring_days = 61:152
    summer_days = 153:245
    fall_days = 246:335
  }
  else{
    wint_days = 1:59
    spring_days = 60:151
    summer_days = 152:244
    fall_days = 245:334
  }
  print("Calculating winter...")
  #
  #
  # Creating seasonal information. May take awhile.
  s_winter1 = life_series_mnth(wint_yr,wint_days_prev,l=l,func=func,norm=norm)
  s_winter2 = life_series_mnth(d,wint_days,l=l,func=func,norm=norm)
  s_winter = (s_winter1 + s_winter2)/2
  print("Calculating spring...")
  s_spring = life_series_mnth(d,spring_days,l=l,func=func,norm=norm)
  print("Calculating summer...")
  s_summer = life_series_mnth(d,summer_days,l=l,func=func,norm=norm)
  print("Calculating fall...")
  s_fall = life_series_mnth(d,fall_days,l=l,func=func,norm=norm)
  fin = list(s_winter,s_spring,s_summer,s_fall) # creating list
  return(fin)
  # returns a list
  # list [[1]] will be winter
  # list [[2]] will be spring
  # list [[3]] will be summer
  # list [[4]] will be fall
}




######################################################################
#
# Looks at the lifetimes over latitudes
# 
#
# I'll be honest I don't know exactly why i made this
# other than the function that plots it (life_plot)
#
######################################################################
life_time = function(d,day,l,nfl=T,norm=F,peaks=F){
  # d = 3-dimensional matrix
  # day = day of the year
  n = dim(d)[2] #the latitudes to use
  omax = max(d) #finding max
  omin = min(d) #finding min
  peakwut = NULL #Initializing peaks
  wut = NULL # Initializing to keep track of 
  nf = NULL # Initializing (n)umber of (f)eatures - nf
  life_avg = NULL
  for (i in 1:n){
    peaklocs = NULL # Initializing peak locations
    m = d[,i,day] # temporary vector
    if (norm == T){ # flag to check if normalizing
      m = (m - omin)/(omax - omin)
    }
    tmp = gridDiag(FUNvalues=fix_series(m),sublevel=FALSE)[['diagram']] # Morse filtration
    b = tmp[,3] # Birth
    lif = tmp[,3] - tmp[,2] # Lifetime
    lif = lif[lif > l] # Lifetime cutoff
    life_avg = c(life_avg,sum(lif)) # Again, it's sum even though named life_avg
    nf = c(nf,length(lif)) # counting number of features
    lif = c(lif,rep(0,20-length(lif))) # creating zeros to have constant vector length
    wut = cbind(wut,lif) # binding together
    b = b[lif > l] # birth that meets the cutoff
    if (length(b) >= 1){ # if we have more than one birth meeting requirement
      for (j in 1:length(b)){
        loc = which(b[j] == m) # Finding the location of the peak relative to the vector
        if (length(loc) == 1){ # checks if there is more than one location for that peak
          peaklocs = c(peaklocs,loc) # sets peak
        }
        else{ #if more than one location for the peak is found
          found = F # setting flag to F
          for (k in 1:length(loc)){ # breaking if we find the peak
            if (found == T){
              break
            }
            if (loc[k] == 1 | loc[k] == dim(d)[1]){ # This checks borders if we land on them
              loc[k] = 2
            }
            if (m[loc[k]] >= m[loc[k]-1] & m[loc[k]] >= m[loc[k]+1]){ # checks if the peak is the "top"
              peaklocs = c(peaklocs,loc[k]) # sets peak
              found = T # set that we found peak
            }
          }
        }
      }
    }
    peaklocs = c(peaklocs,rep(0,20-length(peaklocs))) # fixing format
    peakwut = cbind(peakwut,peaklocs) # putting together
  }
  fin=cbind(nf,life_avg) # creating list
  if (peaks == T){
    return(peakwut)
  }
  if (nfl == T){
    return(fin)
  }
  else{
    return(wut) 
  }
}




######################################################################
#
# This plots the "hill plots" that I use in my thesis
# It shows individual waves and their lifespans and stacks them
# It also overlays information (mean, sd, etc...)
#
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
#
#
######################################################################
life_plot = function(d,day,l,title="",norm=F,overlay=mean,startlat){
  # d = year like y2020, y2017, etc... 3D matrix
  # day = day to plot
  # l = lifetime to ignore
  # title = Title on the plot
  # norm = Whether we normalize or not. Only really changes axis
  # overlay = mean, sd, max, min, var, etc...
  #   overlays a dotted line to look at certain stats
  fc = c(0,"blue","darkorange","cadetblue","deeppink","burlywood1",
         "cornflowerblue","red")
  mn = life_time(d,day,l=l,nfl=F,norm=norm,peaks=F)
  ind = NULL
  for (i in 1:dim(mn)[1]){
    if (sum(mn[i,]) == 0){
      ind = c(ind,i)
    }
  }
  n = ncol(d)
  mn = mn[-ind,]
  fig = plot_ly(y=mn[1,],x=seq(startlat,startlat+(n-1)*2.5,by=2.5),type='scatter',mode='none',
                stackgroup='one',fillcolor='tan',
                name="Feature 1")
  for (i in 2:dim(mn)[1]){
    fig = fig %>% add_trace(y = mn[i,],fillcolor=fc[i],
                            name=paste("Feature",i,sep=" "))
  }
  mn[mn == 0] = NA
  mns = apply(mn,2,overlay,na.rm=TRUE)
  f = deparse(substitute(overlay))
  fig = fig %>% add_trace(y = mns,type='scatter',mode='lines',stackgroup='two',
                          line = list(width=2, dash='dash', color='black'),
                          name = f)
  fig = fig %>% layout(xaxis = list(title = "Latitude N"),
                       yaxis = list(title = "Feature Height"),
                       title = title)
  colnames(mn) = paste(seq(startlat,startlat+(n-1)*2.5,by=2.5),"N",sep=" ")
  print(mn)
  print(colSums(!is.na(mn)))
  return(fig)
}



######################################################################
#
# peak_plot measures where the features are detected and attempts
# to plot them in an animation so we can see the movement of the waves
# it does this on normalized values
# so set l in some range of (0,0.1)
#
#
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
#
######################################################################
peak_plot = function(d,l,startlat){
  n = dim(d)[3]
  what = NULL
  y = rep(1:320,365)
  for (i in 1:n){
    print(i)
    tmp = life_time(d,i,l=l,nfl=F,norm=T,peaks=T)
    tmp = as.vector(tmp)
    what = c(what,tmp)
  }
  what[what == 0] = NA
  days = rep(1:365,1,each=320)
  what = cbind(y,what,days)
  figgy = as.data.frame(what) %>% plot_ly(x=~what,frame=~days,
                                          type = 'scatter', mode = 'markers')
  figgy = figgy %>% animation_opts(125,transition=0,redraw=TRUE)
  figgy = figgy %>% layout(xaxis = list(title = "Longitude"),
                           shapes = list(type = 'rect',
                                         fillcolor='green',opacity = 0.3,
                                         x0 = 112,x1 = 137,
                                         xref = "x",
                                         y0 = 0, y1 = 70,
                                         yref = "y"),
                           yaxis = list(range = c(15,60),
                                        tickvals = seq(15,60,by=2.5),
                                        ticktext=seq(startlat,startlat+(n-1)*2.5,by=2.5),
                                        title = 'Latitude'))
  return(figgy)
}



######################################################################
#
# For this function to work, your 30 year climatology must be named
# y30yr, and be on the same pressure level with the same lat/lon
#
# It can compare various measurements.
# For example, comp_mnth(y2010,l=20,mean,"jan") compares the means in 
#    the month of January in 2010 to the means of the base 1981-2010
#    January climatology. You can specify days if you want that particular
#    day compared
#
#
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
#
#
# this is one of more the more complicated ones
#
######################################################################
comp_mnth = function(d,l,func,mnth,life=T,norm=F,day=NULL){
  #
  # This function compares any particular month OR day to
  # the 30 year climatology
  # d = The year. Formatted yYEAR
  #     For example, 2010 is y2010
  # l = the "cutoff" if we want to ignore certain small lifetimes
  #     if norm=T, this should be between around 0-0.1
  #     if you go higher than that it might throw errors
  # func = function to evaluate the information
  #        mean, median, sd, max, min
  # mnth = Month to look at
  #        accepts string for first three letters of month
  #        For example, March is "mar" and June is "jun"
  # life = T compares lifetimes of features put through the function
  #        if F, it counts features (takes longer)
  # norm = F evaluates information as raw data
  #        T normalizes it over the YEAR. So it uses max(year) and min(year)
  # day = Day of the month to look at. Defaults to NULL if you want to look
  #       at the entire month
  yr = str_remove(deparse(substitute(d)),"y")
  if (norm == T){
    d = (d - min(d))/(max(d) - min(d))
    y30yr_n = (y30yr - min(y30yr))/(max(y30yr) - min(y30yr))
  }
  leap = F
  ndays = dim(d)[3]
  if (ndays == 366){
    leap = T
  }
  days = switch(mnth,
                "jan" = 1:31,
                "feb" = 32:59,
                "mar" = 60:90,
                "apr" = 91:120,
                "may" = 121:151,
                "jun" = 152:181,
                "jul" = 182:212,
                "aug" = 213:243,
                "sep" = 244:273,
                "oct" = 274:304,
                "nov" = 305:334,
                "dec" = 335:365
  )
  if (!is.null(day)){
    days = min(ind) + day - 1
  }
  if (life == T){
    if (norm == T){
      base = life_series_mnth(y30yr_n,days,l,func)
    }
    else{
      base = life_series_mnth(y30yr,days,l,func)
    }
  }
  else if (life == F){
    base = qk_series(y30yr[,,days],l,func)
  }
  if (leap == T & mnth != "jan" & mnth != "feb"){
    ind = ind + 1
  }
  if (leap == T & mnth == "feb"){
    ind = c(ind,60)
  }
  if (life == T){
    if (norm == T){
      comp = life_series_mnth(d,days,l,func)
    }
    else{
      comp = life_series_mnth(d,days,l,func)
    }
  }
  else if (life == F){
    comp = qk_series(d[,,days],l,func)
  }
  if (!is.null(day)){
    base[which(base == 0)] = NA
    comp[which(comp == 0)] = NA
    base = apply(base,2,func,na.rm=TRUE)
    comp = apply(comp,2,func,na.rm=TRUE)
  }
  fig = plot_ly(y=base,type="violin",name=paste(mnth,"climatology",sep="-"),
                box = list(visible = T),
                meanline = list(visible = T),
                color = I("turquoise1"))
  fig = fig %>% add_trace(y=comp,type="violin",name=paste(mnth,yr,sep="-"),
                          box = list(visible = T),
                          meanline = list(visible = T),
                          color = I("springgreen3"))
  f=deparse(substitute(func))
  if (norm == T){
    f = paste("Normalized",f,sep="-")
  }
  if (life == T){
    nm = "life"
  }
  else{
    nm = "features"
  }
  if (!is.null(day)){
    f = paste(f,nm,"on day",day,"of",mnth,sep=" ")
  }
  else{
    f = paste(f,nm,"in",mnth,sep=" ")
  }
  fig = fig %>% layout(yaxis = list(title = f))
  return(fig)
}

######################################################################
#
# This is to just plot the anomaly plot
#
#
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
# THIS FUNCTION NEEDS MORE COMMENTS STILL
#
#
######################################################################

plot_anom = function(d,day){
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366)
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12)
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December")
  mnth = month_labels[findInterval(day,month_intervals)]
  mnth_str = month[findInterval(day,month_intervals)]
  day_num = day - month_intervals[findInterval(day,month_intervals)] + 1
  ref = avg_sub[,,mnth]
  anom = d[,,day] - ref
  anom = rbind(anom[73:144,],anom[1:72,]) #fixing how it looks
  
  rownames(anom) = paste("x",1:144,sep="")
  colnames(anom) = paste("y",1:16,sep="") #naming things
  
  anom = melt(anom) #Meltiiiiiiing. What a world what a world
  long = paste("x",c(1,seq(0,144,by=8)),sep="")
  long = long[-2]
  latg = paste("y",seq(1,16,by=1),sep="")
  longlab = seq(0,144,by=8)*2.5 - 180
  latlab = seq(42.5,80,by=2.5)
  
  p1 = ggplot(anom,aes(x = Var1, y = Var2)) +
    geom_raster(aes(fill = value),interpolate = TRUE) + 
    scale_fill_gradient2(low="blue",mid="white",high="red",
                         midpoint = 0, limits = range(anom$value)) + 
    theme_classic() +
    scale_x_discrete(name = "Longitude",
                     breaks = long, labels = as.character(longlab)) +
    scale_y_discrete(name = "Latitude",
                     breaks = latg, labels = as.character(latlab)) +
    ggtitle(paste("Anomaly plot for",mnth_str,day_num,sep=' '))
  print(p1)
  return(anom)
}

# Makes the anomaly plots for the whole year

make_all_anom = function(d){
  all_days = dim(d)[3]
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366)
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12)
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December")
  for (i in 1:all_days){
    mnth = month_labels[findInterval(i,month_intervals)]
    mnth_str = month[findInterval(i,month_intervals)]
    day_num = i - month_intervals[findInterval(i,month_intervals)] + 1
    png(paste(mnth_str,day_num,".png",sep=""), width = 1057, height = 556)
    plot_anom(d,i)
    dev.off()
  }
}

make_all_seahorse = function(d,startlat){
  all_days = dim(d)[3]
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366)
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12)
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December")
  for (i in 1:all_days){
    print(i)
    mnth = month_labels[findInterval(i,month_intervals)]
    mnth_str = month[findInterval(i,month_intervals)]
    day_num = i - month_intervals[findInterval(i,month_intervals)] + 1
    title = paste("Seahorse for",mnth_str,day_num,sep=" ")
    png(paste("SH_",mnth_str,day_num,".png",sep=""), width = 1057, height = 556)
    count_gd_day(d,i,plot=T,title=title,l=0,startlat=startlat)
    dev.off()
  }
}

make_all_horsiez = function(d,startlat){
  all_days = dim(d)[3]
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366)
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12)
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December")
  for (i in 1:all_days){
    print(i)
    mnth = month_labels[findInterval(i,month_intervals)]
    mnth_str = month[findInterval(i,month_intervals)]
    day_num = i - month_intervals[findInterval(i,month_intervals)] + 1
    title = paste("HorsieZ for",mnth_str,day_num,sep=" ")
    png(paste("hZ_",mnth_str,day_num,".png",sep=""), width = 1057, height = 556)
    horsieZ(d,daycomp = i,mtx = america_anom30[[mnth]],plotting=T,title=title,l=0,startlat=startlat)
    dev.off()
  }
}


######################################################################
#
# This is to make a series of plots
# And will likely be more annoying to look at
#
# I was right. The colorings are very annoying to get correct
# I need to create a list of mins/maxs for each levels
#
######################################################################

# These are all the levels we're possibly interested in
load(file = "north30_1000mb.RData")
load(file = "north30_850mb.RData")
load(file = "north30_700mb.RData")
load(file = "north30_500mb.RData")
load(file = "north30_300mb.RData")
load(file = "north30_200mb.RData")
load(file = "north30_100mb.RData")
load(file = "mm_yrs.RData") # this is a matrix of all mins/maxes

plot_levels = function(year,day,lats,startlat){
  #
  # Lats should be the exact same as when you created the year
  # EX: if you made y2021 = s2021[,1:13,,] then lats should be 1:13
  # startlat should be the actual measurement of your starting latitude
  # if you start at 0, should be 0
  # if you start at 42.5N, should be 42.5
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366) # setting intervals
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12) # setting month labels
  mnth = month_labels[findInterval(day,month_intervals)] # finding the correct month
  d.yr = nc_open(paste("hgt",year,"nc",sep='.')) %>%
    ncvar_get("hgt") # fetching the correct year to compare
  d.yr = d.yr[,lats,,] # setting the correct latitudes
  par(mar = c(2.5, 0.5, 3.5, 0.5)) # margins -- THIS CAN BE EDITED
  par(mfrow = c(4,2)) # putting plots on the same place
  
  #
  #
  #
  #
  # This creates all the level plots
  horsieZ(d.yr[,,1,],day,mtx=north30_1000mb[[mnth]],plotting=T,startlat=startlat,l=0,title="1000mb",lvls=mm_yrs[1,])
  horsieZ(d.yr[,,3,],day,mtx=north30_850mb[[mnth]],plotting=T,startlat=startlat,l=0,title="850mb",lvls=mm_yrs[2,])
  horsieZ(d.yr[,,4,],day,mtx=north30_700mb[[mnth]],plotting=T,startlat=startlat,l=0,title="700mb",lvls=mm_yrs[3,])
  horsieZ(d.yr[,,6,],day,mtx=north30_500mb[[mnth]],plotting=T,startlat=startlat,l=0,title="500mb",lvls=mm_yrs[4,])
  horsieZ(d.yr[,,8,],day,mtx=north30_300mb[[mnth]],plotting=T,startlat=startlat,l=0,title="300mb",lvls=mm_yrs[5,])
  horsieZ(d.yr[,,10,],day,mtx=north30_200mb[[mnth]],plotting=T,startlat=startlat,l=0,title="200mb",lvls=mm_yrs[6,])
  horsieZ(d.yr[,,12,],day,mtx=north30_100mb[[mnth]],plotting=T,startlat=startlat,l=0,title="100mb",lvls=mm_yrs[7,])
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December") # Labels for months
  day_num = day - month_intervals[findInterval(day,month_intervals)] + 1 #day of the month
  title(paste(month[mnth],day_num,year),line=-2,outer=TRUE) #creating a title
}

plot_levels_year = function(year,lats,startlat){ # This does the same thing as above but for a range of days
  #
  #
  # saves an image with the format lvls-Month-DayOfMonth
  # For example, april 1st is lvlsApril1.png
  d.yr = nc_open(paste("hgt",year,"nc",sep='.')) %>%
    ncvar_get("hgt") # fetching the correct year
  d.yr = d.yr[,lats,,] # temporary
  n = dim(d.yr)[4] # number of days
  month_intervals = c(1,32,60,91,121,152,182,213,244,274,305,335,366) #months
  month_labels = c(1,2,3,4,5,6,7,8,9,10,11,12) #labels
  month = c("January","February","March","April",
            "May","June","July","August",
            "September","October","November","December") # string labels
  for (i in 1:n){ #loop to look over all days
    mnth = month_labels[findInterval(i,month_intervals)] # another month
    mnth_str = month[findInterval(i,month_intervals)] # month string chosen from month
    day_num = i - month_intervals[findInterval(i,month_intervals)] + 1 #day
    png(paste("lvls",mnth_str,day_num,".png",sep=''),width=1920,height=1080) #saving the image
    plot_levels(year,i,lats=lats,startlat=startlat) #plotting from above
    dev.off() #boop
  }
}
