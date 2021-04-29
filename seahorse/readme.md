This README will be updated in short bursts to explain the individual R scripts that I've added here and what each of them do. There's a fair amount of random things I did in my exploration so this may take a little bit of time. As of right now, there are comments in each of them that gives a short explanation as to what they do.

# seahorse.R

*in this folder, there is a file named bymonth30.RData. This is a 30 year base built with build_comp with l=50. It can be loaded in R with the command load(file="bymonth30.RData"), which will create a variable named bymonth30 that contains all 12 months. To access a month, use bymonth30[[i]] where i is the number of the month, so bymonth30[[2]] would be February. In addition, each month contains 3 subsets of data. Mean is a cell-wise mean of all samples, SD is a cell-wise standard deviation of all samples, and Raw is the raw data used to compute these*

**fix_series(d)**: This function moves the minimums of the curve to the borders so that our Morse filtration does not count features on the border twice, since the data we are working with is technically a loop that goes around the globe. The input is a single 1-D curve (**d**).

**count_gd(d,ln,l)**: Uses Morse filtration at a particular filtration height. This is the basis of the seahorse / horsieZ. It works on 1-D curves and look at only one filtration height. It takes inputs **d**(1-D curve), **ln**(a single filtration height), and **l**(the smallest lifetime to be considered a feature).

**count_gd_heights(d,day,lts,l,omax,nheight)**: Performs count_gd over a range of heights that we can define the "fineness" of. This is the function that creates each "slice" of latitude that we display on the seahorse plot. It takes inputs **d**(3-D Matrix), **day**(day of the year), **lts**(Latitudes to look at), **l**(the smallest lifetime to be considered a feature), **omax**(this is mostly to set boundaries for horsieZ and not used otherwise), and **nheight**(the number of "slices" to split filtration heights into).

**count_gd_day(d,day,plot,title,l,summary,omax,nheight)**: Performs count_gd_day over all possible latitudes for the matrix given and plots the seahorse. It takes inputs **d**(3-D Matrix), **day**(day of the year), **plot**(T/F if we should plot the output or not), **title**(title on the plot), **l**(the smallest lifetime to be considered a feature), **summary**(T/F to display some overlaid lines of means and such, kind of messy and may be removed), **omax**(this is mostly to set boundaries for horsieZ and not used otherwise), and **nheight**(the number of "slices" to split filtration heights into).

**build_comp(daycomp,year_comp,days_back,years_back,l)**: Builds a baseline to compare against. Each pixel will have (days_back+1) by years_back number of samples. For example, days_back=30 and years_back=30 would contain 930 samples per pixel. If your matrix is 100by100, then you would have 10,000 pixels, each with 930 samples, so this is a relatively large process for less powerful laptops like mine :(. It takes inputs **daycomp**(the day we "start at" in the year, 1-365), **year_comp**(the year we "start at", ex 2010 or 2015), **days_back**(the number of days to step backwards from daycomp. If daycomp = 31 and days_back=30, we do days 31-1), **years_back**(the number of years to step backwards from year_comp),**l**(the smallest lifetime to be considered a feature)

**horsieZ(d,daycomp,mtx,l,plotting,title)**: This function compares a particular day to the base that is built with build_comp. It compares them via computing a Z-score for each pixel and displays them in a similar fashion to the seahorse plot from count_gd_day. Hopefully, this helps identify extreme weather patterns. It takes input **d**(data for the year we are wanting to use as a comparison. format of y2010, y2015, etc...), **daycomp**(the day we are going to compare), **mtx**(the base we created in build_comp), **l**(the smallest lifetime to be considered a feature), **plotting**(T/F if we want to display the plot, **title**(title of the plot).



# seahorse_misc.R

**will add information later**

# horsieZ_misc.R

**will add information later**
