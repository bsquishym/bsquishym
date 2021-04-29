# Requires having all the information in seahorse.R active

# These functions are mostly small explorations into the horseiZ function
# It is currently unknown if any of this is useful, but there is no harm
# in curiousity according to my cat



######################################################################
#
# This little fella goes through each horsieZ for each day that you list
#   in days. It then compares each entry to our extreme chebyshev z of 3.
#   it outputs a vector listing how many entries were considered extreme.
#   It displays a series of these counts
#
######################################################################
count_rej_days = function(d,days,base,l,title="",lats=1:16){
  #
  # d = yearly data. ex y2010, y2020, in that format
  # days = a range of days in 1-365
  # base = base from build_comp in seahorse.R
  # l = lifetime cutoff
  # lats = if you want to look at specific latititudes. Defaults to all of ours
  #
  # For example, count_rej_days(y2020,days=1:31,base=bymonth30[[1]],l=50) would
  #    create a horsieZ for each day of January in 2020, count the number of
  #    extreme features detected, record that number, and display a series
  #    with hopes of being able to identify particular extreme events
  crit99 = 3
  count99 = NULL
  count95 = NULL
  for (i in days){
    check = horsieZ(d,daycomp=i,base,l=l,plotting=F)
    count99 = c(count99,length(which(abs(check[lats,]) > crit99)))
  }
  fig = plot_ly(y = count99, type = 'scatter', mode = 'lines',name="Extreme Observations")
  fig = fig %>% layout(title = title, showlegend = T,
                       xaxis = list(ticktext = days,
                                    title = "Days"),
                       yaxis = list(title = "Rejected Features"))
  print(fig)
  return(count99)
}



######################################################################
#
# This is the big brother of count_rej_days
# This not only displays rejected counts per day but separates it
#   by latitudes as well
# This is displayed on a counter plot via plot_ly
#
######################################################################
count_rej_sh = function(d,days,base,l,title=""){
  #
  # d = yearly data. ex y2010, y2020, in that format
  # days = a range of days in 1-365
  # base = base from build_comp in seahorse.R
  # l = lifetime cutoff
  #
  # For example, count_rej_sh(y2020,days=1:31,base=bymonth30[[1]],l=50) would
  #    create a horsieZ for each day of January in 2020, count the number of
  #    extreme features detected at each latitude, record that number, and display
  #    a contour plot with blobs that show which latitudes the extreme features
  #    showed up at with hopes of being able to identify particular
  #    extreme events
  crit99 = 3
  count_tot = NULL
  for (i in days){
    print(paste("Day",i,sep=" "))
    count99 = NULL
    check = horsieZ(d,daycomp=i,base,l=l,plotting=F)
    for (j in 1:nrow(check)){
      count99 = c(count99,length(which(abs(check[j,]) > crit99)))
    }
    count_tot = cbind(count_tot,count99)
  }
  figgy = plot_ly(z = count_tot, type='contour',name='Rejected Features')
  figgy = figgy %>% layout(xaxis = list(title = "Day"),
                           yaxis = list(tickvals = seq(0,16,by=1),
                                        ticktext=seq(42.5,80,by=2.5),
                                        title = 'Latitude'),
                           title = title)
  print(figgy)
  fin = NULL
  colnames(count_tot) = paste("Day",seq(min(days),max(days),1),sep="")
  return(count_tot)
}




######################################################################
#
# This is the same as count_rej_sh except the y-axis counts rejections
#    based off of geopotential height with the hopes that we can see
#    rejected features in geopotential heights that are usually not active
#
#
######################################################################
count_rej_shGPH = function(d,days,base,l,title=""){
  #
  # d = yearly data. ex y2010, y2020, in that format
  # days = a range of days in 1-365
  # base = base from build_comp in seahorse.R
  # l = lifetime cutoff
  #
  # For example, count_rej_shGPH(y2020,days=1:31,base=bymonth30[[1]],l=50) would
  #    create a horsieZ for each day of January in 2020, count the number of
  #    extreme features detected at each geopotential height, record that number, and display
  #    a contour plot with blobs that show which latitudes the extreme features
  #    showed up at with hopes of being able to identify particular
  #    extreme events
  crit99 = 3 #chebyshev extreme value
  count95 = NULL
  count_tot = NULL
  for (i in days){
    count99 = NULL
    check = horsieZ(d,daycomp=i,base,l=l,plotting=F) #creating horsieZ
    for (j in 1:ncol(check)){
      count99 = c(count99,length(which(abs(check[,j]) > crit99))) #vs chebyshev
    }
    count_tot = cbind(count_tot,count99) #putting together
  }
  yax = seq(round(4594+l,-2),round(6018+l,-2),by=100) #y-axis
  figgy = plot_ly(z = count_tot, type='contour',name='Rejected Features')
  figgy = figgy %>% layout(xaxis = list(title = "Day"),
                           yaxis = list(tickvals = seq(0,100,length.out=16),
                                        ticktext=yax,
                                        title = 'Geopotential Height'),
                           title = title)
  print(figgy)
  fin = NULL
  return(count_tot)
}




######################################################################
#
# This set of 3 functions is for looking at rarely occuring features
# tally_zero and tally_horse just count
# pls excuse my variable names they're pretty strange
# The main one is rare_track
#
######################################################################
tally_zero = function(v){
  k = 0
  for (i in 1:length(v)){
    if (v[i] != 0){
      k = k + 1
    }
  }
  return(k)
}

tally_horse = function(d){
  well = apply(d,c(1,2),tally_zero)
}

rare_track = function(d){
  #
  # d in this case is a raw version of one of the bases
  # for example, d = bymonth30[[1]]$Raw, the raw data of Januaries
  #
  # the function goes through the raw data by year and looks for
  #   times that a feature occured where it never had occured before
  #   it then inputs that year into a tracking matrix to display the
  #   year the it was first observed
  #
  # this could be expanded to a larger base for a more thorough look
  # like making a 40 year horsieZ base and tracking throughout
  # since on the 30 year base, rare events still occur up until 2010
  # which is the last year we look
  #
  # returns a matrix with entries 1981-2010 which is the year activity
  #   was first observed
  n = dim(d)[3]
  track.matrix = matrix(nrow=16,ncol=100,0) #tracking matrix
  what = seq(n,0,by=-31) #sequence of years to check
  what[length(what)] = 1
  for (i in 1:(length(what)-1)){
    look = seq(what[i],what[i+1],by=-1) #looking through years
    wee = tally_horse(d[,,look])
    ind = which(wee == 1) #checking rare events
    for (k in ind){
      if (track.matrix[k] == 0){ #comparing rare events to see if it's the first one
        track.matrix[k] = 1980 + i #setting to year if it's the first time
      }
    }
  }
  return(track.matrix) #matrix of year entries
}
