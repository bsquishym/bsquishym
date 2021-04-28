This R file contains most of what is needed to do comparisons between base data and a particular day.

It contains the following functions:

fix_series -> Helps prevent double counting for features.

count_gd -> Does a single morse filtration at a particular filtration height and counts features above it
         -> Input is a vector
         
count_gd_heights -> performs count_gd over a range of heights
         -> Input is a vector
         
count_gd_day -> performs count_gd_heights over an entire matrix
             -> Input is a matrix
             
build_comp -> Creates a baseline over the desired number of days/years

horsieZ -> Compares a day vs a baseline
