This README will be updated in short bursts to explain the individual R scripts that I've added here and what each of them do. There's a fair amount of random things I did in my exploration so this may take a little bit of time. As of right now, there are comments in each of them that gives a short explanation as to what they do.

# seahorse.R

*in this folder, there is a file named bymonth30.RData. This is a 30 year base built with build_comp with l=50. It can be loaded in R with the command load(file="bymonth30.RData"), which will create a variable named bymonth30 that contains all 12 months. To access a month, use bymonth30[[i]] where i is the number of the month, so bymonth30[[2]] would be February. In addition, each month contains 3 subsets of data. Mean is a cell-wise mean of all samples, SD is a cell-wise standard deviation of all samples, and Raw is the raw data used to compute these*

*there's also bymonth30_0.RData which is a base with l=0*

**fix_series(d)**: This function moves the minimums of the curve to the borders so that our Morse filtration does not count features on the border twice, since the data we are working with is technically a loop that goes around the globe. The input is a single 1-D curve (**d**).

**count_gd(d,ln,l)**: Uses Morse filtration at a particular filtration height. This is the basis of the seahorse / horsieZ. It works on 1-D curves and look at only one filtration height. It takes inputs **d**(1-D curve), **ln**(a single filtration height), and **l**(the smallest lifetime to be considered a feature).

**count_gd_heights(d,day,lts,l,omax,nheight)**: Performs count_gd over a range of heights that we can define the "fineness" of. This is the function that creates each "slice" of latitude that we display on the seahorse plot. It takes inputs **d**(3-D Matrix), **day**(day of the year), **lts**(Latitudes to look at), **l**(the smallest lifetime to be considered a feature), **omax**(this is mostly to set boundaries for horsieZ and not used otherwise), and **nheight**(the number of "slices" to split filtration heights into).

**count_gd_day(d,day,plot,title,l,summary,omax,nheight)**: Performs count_gd_day over all possible latitudes for the matrix given and plots the seahorse. It takes inputs **d**(3-D Matrix), **day**(day of the year), **plot**(T/F if we should plot the output or not), **title**(title on the plot), **l**(the smallest lifetime to be considered a feature), **summary**(T/F to display some overlaid lines of means and such, kind of messy and may be removed), **omax**(this is mostly to set boundaries for horsieZ and not used otherwise), and **nheight**(the number of "slices" to split filtration heights into).

**build_comp(daycomp,year_comp,days_back,years_back,l)**: Builds a baseline to compare against. Each pixel will have (days_back+1) by years_back number of samples. For example, days_back=30 and years_back=30 would contain 930 samples per pixel. If your matrix is 100by100, then you would have 10,000 pixels, each with 930 samples, so this is a relatively large process for less powerful laptops like mine :(. It takes inputs **daycomp**(the day we "start at" in the year, 1-365), **year_comp**(the year we "start at", ex 2010 or 2015), **days_back**(the number of days to step backwards from daycomp. If daycomp = 31 and days_back=30, we do days 31-1), **years_back**(the number of years to step backwards from year_comp),**l**(the smallest lifetime to be considered a feature)

**horsieZ(d,daycomp,mtx,l,plotting,title)**: This function compares a particular day to the base that is built with build_comp. It compares them via computing a Z-score for each pixel and displays them in a similar fashion to the seahorse plot from count_gd_day. Hopefully, this helps identify extreme weather patterns. It takes input **d**(data for the year we are wanting to use as a comparison. format of y2010, y2015, etc...), **daycomp**(the day we are going to compare), **mtx**(the base we created in build_comp), **l**(the smallest lifetime to be considered a feature), **plotting**(T/F if we want to display the plot, **title**(title of the plot).



# seahorse_misc.R

**qk_series(d,l,func)**: This function creates a series for any particular range of days. The series is a a summary statistic over a seahorse plot created over each day. For example, qk_series(y2020[,,1:31],l=0,mean) will create a series of all days of January, where each entry is the mean of the feature counts of that day. Takes input **d**(year and day in the form of yXXXX[,,range of days]), **l**(the smallest lifetime to be considered a feature), **func**(a function to run over all entries of the seahorse plot, like mean, max, sd, etc...).

**count_gd_season(d,l,func,norm)**: This function looks at a year and subsets it into seasons. If we look at 2010, winter includes December of 2009. We define the seasons as winter (December - February), spring (March - May), summer (June - August), and fall (September - November). It returns a list with 4 entries with [1] being winter and [4] being fall. It takes input **d**(the "reference" year), **l**(the smallest lifetime to be considered a feature), **func**(a function to run over all entries of the seahorse plot, like mean, max, sd, etc...), and **norm**(whether we normalize before performing the function).

**life_time(d,day,l,norm)**: This function looks at a particular day, performs a Morse filtration over all latitudes, and sums up all the lifetimes of all features, which it returns as a vector of values. It takes inputs **d**(the reference year), **day**(the day we look at, 1-365), **l**(the smallest lifetime to be considered a feature), and **norm**(whether we normalize before performing the function).

**life_series_mnth(d,days,l,func,norm)**: This performs life_time over multiple days and then performs a function element-wise throughout the days. It takes inputs **d**(the reference year), **day**(the day we look at, 1-365), **l**(the smallest lifetime to be considered a feature), **func**(a function to run over all entries of the seahorse plot, like mean, max, sd, etc...), and **norm**(whether we normalize before performing the function).

**life_season(d,l,func,norm)**: This does the same thing count_gd_season does, except with feature lifetimes. It takes input **d**(the "reference" year), **l**(the smallest lifetime to be considered a feature), **func**(a function to run over all entries of the seahorse plot, like mean, max, sd, etc...), and **norm**(whether we normalize before performing the function).

**life_time(d,day,l,nfl,norm,peaks)**: This function is a *helper* function for life_plot. It has a few different options. If **nfl = T**, it returns a matrix that contains the number of features(first column) and the sum of their lifetimes. If **nfl = F**, it returns a matrix that contains the lifetimes of each feature in an array. If **peaks = T**, it returns a similar array as nfl but the location of the peaks detected in the curve. It is mostly experimental. The rest of the inputs are the usual **d**(the "reference" year), **day**(the day we look at, 1-365), **l**(the smallest lifetime to be considered a feature), and **norm**(whether we normalize before performing the function).

**life_plot(d,day,l,title,norm,overlay)**: This function displays the wave features we detect with life_time. It stacks the features in an attempt to visualize their relative size and how much of a wave is at a particular latitude. It takes input **d**(the "reference" year), **day**(the day we look at, 1-365), **l**(the smallest lifetime to be considered a feature), **title**(title of the plot), **norm**(whether we normalize before performing the function), and **overlay**(what summary statistic to overlay on the plot as a dotted line).

**peak_plot(d,l)**: This function attempts to plot the peaks on a moving plot to investigate how the waves move. It piggybacks off of the life_time function. It takes only two inputs, **d**(the reference year) and **l**(the smallest lifetime to be considered a feature).

**comp_mnth(d,l,func,mnth,life,norm,day)**: This function directly compares a month or day to the 30-year climatology from NOAA's database. It can compare features or lifetimes. It takes inputs **d**(the reference year), **l**(the smallest lifetime to be considered a feature), **func**(a function to run over all entries of the seahorse plot, like mean, max, sd, etc...), **mnth**(month to use as comparison, first 3 letters of a month for example "jun" for june), **life**(T/F if the plot should be of features or lifetimes), **norm**(whether we normalize before performing the function), and **day**(day to compare. if left NULL, it compares the month).


# horsieZ_misc.R

**will add information later**
