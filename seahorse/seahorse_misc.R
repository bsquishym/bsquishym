# Requires having all the information in seahorse.R active

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
  #   mean will take the mean (over the non-zero feature counts)
  #   median will take the median over the non-zero feature counts
  #   etc...
  srs = NULL
  if (is.na(dim(d)[3])){
    d = abind(d,d,along=3)
    srs = count_gd_day(d,1,plot=F,l=0)
  }
  else{
    n = dim(d)[3] #the days of the 3-dimensional dataset
    for (i in 1:n){
      tmp = count_gd_day(d,i,plot=F,l=0)
      srs = c(srs,func(tmp))
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
  ref_yr = str_remove(deparse(substitute(d)),"y")
  wint_yr_ind = as.numeric(ref_yr) - 1
  wint_yr = eval(parse(text = paste("y",wint_yr_ind,sep="")))
  if (dim(wint_yr)[3] == 366){
    wint_yr = wint_yr[,,336:366]
    d_wint = d[,,1:59]
    d_spring = d[,,60:151]
    d_summer = d[,,152:244]
    d_fall = d[,,245:334]
  }
  else{
    wint_yr = wint_yr[,,335:365]
  }
  if (norm == T){
    dmax = max(d)
    dmin = min(d)
    d = (d - dmin)/(dmax - dmin)
    wint_yr = (wint_yr - dmin)/(dmax - dmin)
  }
  if (dim(d)[3] == 366){
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
  d_winter = abind(wint_yr,d_wint,along=3)
  print("Calculating winter...")
  s_winter = qk_series(d_winter,l=0,func=func)
  print("Calculating spring...")
  s_spring = qk_series(d_spring,l=0,func=func)
  print("Calculating summer...")
  s_summer = qk_series(d_summer,l=0,func=func)
  print("Calculating fall...")
  s_fall = qk_series(d_fall,l=0,func=func)
  fin = list(s_winter,s_spring,s_summer,s_fall)
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
  n = dim(d)[2] #the latitudes to use
  life_avg = NULL
  for (i in 1:n){
    m = d[,i,day]
    if (norm == T){
      omax = max(d)
      omin = min(d)
      m = (m - omin)/(omax - omin)
    }
    tmp = gridDiag(FUNvalues=fix_series(m),sublevel=FALSE)[['diagram']]
    lif = tmp[,3] - tmp[,2]
    lif = lif[lif > l]
    life_avg = c(life_avg,sum(lif))
  }
  return(life_avg)
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
  lives = NULL
  if (is.na(dim(d)[3])){
    d = abind(d,d,along=3)
    lives = life_time(d,1,l)
  }
  else{
    n = dim(d)[3] #number of days in the year
    for (i in days){
      # print(i)
      # print(dim(d))
      # print(deparse(substitute(d)))
      lives = cbind(lives,life_time(d,i,l,norm=norm))
    }
    # print(lives)
  }
  if (is.null(dim(lives))){
    fin = lives
  }
  else{
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
  seasons = c("Winter","Spring","Summer","Fall")
  ref_yr = str_remove(deparse(substitute(d)),"y")
  wint_yr_ind = as.numeric(ref_yr) - 1
  wint_yr = eval(parse(text = paste("y",wint_yr_ind,sep="")))
  if (dim(wint_yr)[3] == 366){
    wint_days_prev = 336:366
    wint_days = 1:59
    spring_days = 60:151
    summer_days = 152:244
    fall_days = 245:334
  }
  else{
    wint_days_prev = 335:365
  }
  if (dim(d)[3] == 366){
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
  print(wint_days_prev)
  s_winter1 = life_series_mnth(wint_yr,wint_days_prev,l=l,func=func,norm=norm)
  s_winter2 = life_series_mnth(d,wint_days,l=l,func=func,norm=norm)
  s_winter = (s_winter1 + s_winter2)/2
  print("Calculating spring...")
  s_spring = life_series_mnth(d,spring_days,l=l,func=func,norm=norm)
  print("Calculating summer...")
  s_summer = life_series_mnth(d,summer_days,l=l,func=func,norm=norm)
  print("Calculating fall...")
  s_fall = life_series_mnth(d,fall_days,l=l,func=func,norm=norm)
  fin = list(s_winter,s_spring,s_summer,s_fall)
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
  omax = max(d)
  omin = min(d)
  peakwut = NULL
  wut = NULL
  nf = NULL
  life_avg = NULL
  for (i in 1:n){
    peaklocs = NULL
    m = d[,i,day]
    if (norm == T){
      m = (m - omin)/(omax - omin)
    }
    tmp = gridDiag(FUNvalues=fix_series(m),sublevel=FALSE)[['diagram']]
    b = tmp[,3]
    lif = tmp[,3] - tmp[,2]
    lif = lif[lif > l]
    life_avg = c(life_avg,sum(lif))
    nf = c(nf,length(lif))
    lif = c(lif,rep(0,10-length(lif)))
    wut = cbind(wut,lif)
    b = b[lif > l]
    if (i == 16){
      print(length(b))
      print(b == 0)
      print(b)
    }
    if (length(b) >= 1){
      for (j in 1:length(b)){
        loc = which(b[j] == m)
        if (length(loc) == 1){
          peaklocs = c(peaklocs,loc)
        }
        else{
          found = F
          for (k in 1:length(loc)){
            if (found == T){
              break
            }
            if (loc[k] == 1 | loc[k] == 144){
              loc[k] = 2
            }
            if (m[loc[k]] >= m[loc[k]-1] & m[loc[k]] >= m[loc[k]+1]){
              peaklocs = c(peaklocs,loc[k])
              found = T
            }
          }
        }
      }
    }
    peaklocs = c(peaklocs,rep(0,20-length(peaklocs)))
    peakwut = cbind(peakwut,peaklocs)
  }
  fin=cbind(nf,life_avg)
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
#
######################################################################
life_plot = function(d,day,l,title="",norm=F,overlay=mean){
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
  mn = mn[-ind,]
  fig = plot_ly(y=mn[1,],x=seq(42.5,80,by=2.5),type='scatter',mode='none',
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
  colnames(mn) = paste(seq(42.5,80,by=2.5),"N",sep=" ")
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
######################################################################
peak_plot = function(d,l){
  n = dim(d)[3]
  what = NULL
  y = rep(1:320,365)
  for (i in 1:n){
    print(i)
    tmp = life_time2(d,i,l=l,nfl=F,norm=T,peaks=T)
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
                           yaxis = list(range = c(0,70),
                                        tickvals = seq(0,70,length.out=16),
                                        ticktext=seq(42.5,80,by=2.5),
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
    base = qk_series(y30yr[,,ind],l,func)
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
    comp = qk_series(d[,,ind],l,func)
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
