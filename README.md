These codes are a work in progress. Two dependencies must be installed 

```
ssc install moremata; ssc install bspline
copy "https://raw.githubusercontent.com/paulbousquet/lpdecomp/main/lpdecomp.ado" lpdecomp.ado
```

The code will create the weights needed to plot the evolution of the weighted average throughout time, as well as plot the IRF. 

# Example Implementation

Using data from [Ramey (2016)](https://econweb.ucsd.edu/~vramey/research/Shocks_HOM_Ramey_published_corrected.pdf)

I also use a custom program (`standshock.ado` in this repo) to standardize the monetary policy shock so the IRFs don't need to be rescaled.

```
drop _all
clear all

import excel Monetarydat.xlsx, sheet("Monthly") firstrow case(l)

gen mdate = m(1959m1) + _n-1
tsset mdate, m
gen t = _n

// standardize shock 
replace ff1 = 0 if abs(ff1) < .005
standshock ff1 


// create controls 
gen zlb = (tffr < .2)

local set lcpi lip ebp 

foreach var of local set {
	gen d`var' = 100 * D.`var'
	local mount `mount' d`var'
}

gen dunemp = D.unemp

local mount dunemp dlip debp

drop if mdate > m(2020m1)

// Run decomposition with following options: 
// Horizon of 2 years
// A years worth of lags
// Don't include lagged zlb for collinearity
// Cumulative IRF
// Don't scale IRF 
lpdecomp dlcpi ff1_std `mount', h(24) lag(11) contemp(zlb) cum noadj
```
This code creates a matrix stored as `makes`, which you can turn into a variable using `svmat double makes, names(makes)`. For more information on the options associated with this function, see my [smooth local projections repo](https://github.com/paulbousquet/SmoothLP). 

Visualization is especially a work in progress. Here's some code to visualize the 0 horizon point estimate. 
```
svmat double makes, names(makes)

// Find where makes1 ends
summarize time if !missing(makes1)
local max_t = r(max)

// Create a shifted time variable for the IRF
gen time_shifted = time + `max_t' + 1

* Create rescaled date for the IRF section
* Create date variable starting from November 1988
gen date = tm(1988m11) + time - 1 if !missing(time)
format date %tm

* Create rescaled date for the IRF section
gen date_rescaled = tm(1988m11) + time_shifted + (time_shifted - time_shifted[1]) * 9 - 1 if time <= 20
format date_rescaled %tm

* Find the last date in the main series to place the vertical line
sum date if !missing(makes1)
local last_main = r(max)
local irf_line = `last_main' + 1

* Get y-coordinate for text placement
sum makes1 result1
local text_y = abs(r(max)) * 3.4  // Place text at 90% of max height

* Plot with extended dates and IRF indicator
twoway (line makes1 date if !missing(makes1), lcolor(blue) lwidth(medthick)) ///
       (rarea irc1 irc2 date_rescaled if time<=20, fcolor(purple%15) lcolor(gs13) lw(none) lpattern(solid)) ///
       (scatter result1 date_rescaled if time<=20, c(l) clp(l) ms(i) clc(black) mc(black) clw(medthick)), ///
       xline(`irf_line', lcolor(black) lwidth(medium)) ///
       text(`text_y' `irf_line' "← IRF begins", place(e) size(small) color(black)) ///
       xlabel(`=tm(1990m1)' `=tm(1995m1)' `=tm(2000m1)' `=tm(2005m1)' `=tm(2010m1)' `=tm(2015m1)' , format(%tmCY)) ///
       xtitle("") ///
	   legend(order(1 "Decomposition" 3 "IRF") position(6) ring(0) cols(1)) ///
	   title("CPI Impulse Response + Decomposition of Period 0 Coefficient") ///
	    yscale(range(-0.2 0.1)) ylabel(-0.2(0.05)0.1)
 

```

To add another line, simply do something like 

```
(line makes7 date if !missing(makes7), lcolor(blue) lwidth(medthick)) ///
```



