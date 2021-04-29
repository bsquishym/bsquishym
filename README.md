Hi. I just recently graduated with my MS in statistics. I figured I needed a place to house all of my random functions and explorations that I encountered in the process of creating my thesis.

All of my files involve atmospheric geopotential height data that can be found at https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.pressure.html which includes geopotential heights for for multiple pressure levels (1000 millibar - 10 millibar).

My thesis involved looking at the structure of atmospheric pressure waves at certain pressure heights (500mb) to try to identify extreme atmospheric activity. It uses Topological Data Analysis to find "features" and compares them vs a baseline that we create using the same process. The hope is that these detected features are able to detect a certain type of wave known as a Rossby wave that is associated with extreme weather. The waves are defined via a spatial frequency (ex. 5-wave has 5 waves radiating downwards from the north pole) and the larger number waves (5-7) are the ones that we are mostly curious about.

Each year's data consists of a 3-dimensional matrix with dimensions (144,73,365) on normal years and (144,73,366) on leap years. We have 144 longitudes, starting at 0deg - 357.5deg at 2.5deg intervals (which translates to (0,180e,180w,0) when we write it in terms of E/W, and 73 latitudes from 90deg S - 90deg N. The dimension with 365 matrix slices is the day of the year. We concentrate on the latitudes 42.5N - 80N.

We take each day, fix a particular latitude, and run a perform a Morse filtration over the 1-d curve. Morse filtration can be thought of as "draining" a bathtub and keeping up when curves are above the waterlevel. We choose a range of filtration heights to perform this process, and record the number of features above each level. This process is possible over the entire matrix, but due to the nature of the filtration, we lose some location data. Therefore, we fix the latitudes as we expect the Rossby waves to manifest longitudinally. We do this for each latitude in our band and display the number of features on a color coded plot that we call a seahorse plot.

We then create a baseline to compare this information to. We take the years 1981-2010, which is the norm for calculating "normal" weather, and perform this morse filtration over each day. For example, we take every January from 1981-2010 and create a seahorse plot. Each filtration height-latitude combination contains 930 samples of the weather activity at that height-latitude. We take those 930 samples and calculate the mean and standard deviation. We then take a day we want to compare, say January 1st, 2020, and create the seahorse plot for that particular day. Then, we compute the Z-score for each "pixel." Since the distributions of these feature samples is not normal, we interpret them via Chebyshev's rule and consider Z-scores of |0-2| to be normal, |2-3| to be unusual, and |3+| to be extreme. This allows us to give a numerical representation of the atmospheric activity.

While this was the subject of my thesis, I am still in collaboration with my advisor Dr. Lynne Seymour and Dr. Adam Jaeger as this project has become very interesting to me and I'd like to see how the project evolves, even if I am not required to be involved in it.

<!---
bsquishym/bsquishym is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
