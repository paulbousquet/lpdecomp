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

import excel using ///
"https://raw.githubusercontent.com/paulbousquet/lpdecomp/main/Monetarydat.xlsx", ///
sheet("Monthly") firstrow case(lower)

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
// Don't include lagged zlb (collinearity)
// Cumulative IRF
// Don't scale IRF 
lpdecomp dlcpi ff1_std `mount', h(24) lag(11) contemp(zlb) cum noadj decompgraph
```
By default, this will show the decomposition for the impact horizon. You can choose a different horizon by adding the option `makescol(k)`. 

Note this code assumes your `tsset` in Stata is done with respect to an actual time variable. 

For more customizability, this code creates a variable using `svmat double makes, names(makes)`. For more information on the options associated with this function, see my [smooth local projections repo](https://github.com/paulbousquet/SmoothLP). 

Here are the relevant parts of the function to recreate graphs exactly which you can then customize further 
```
qui sum makes`makescol' irc2
	local y_min = r(min)
	local y_max = r(max)
	local y_range = `y_max' - `y_min'
	local text_y = (`y_max' - result1[1])/2 
	
	// Build xlabel values: 4-year increments starting from even year
	local first_year = year(dofm(`first_main'))
	if mod(`first_year', 4) != 0 {
		local first_year = `first_year' + (4 - mod(`first_year', 4))
	}
	local last_year = year(dofm(`last_main'))
	local xlabs ""
	forvalues yr = `first_year'(4)`last_year' {
		local lab_date = tm(`yr'm1)
		local xlabs "`xlabs' `lab_date'"
	}
	
	// Text position slightly to the right of the line
	local text_x = `irf_line' + 5

twoway (line makes`makescol' makes_date if !missing(makes`makescol'), lcolor(blue) lwidth(medthick)) ///
				   (rarea irc1 irc2 date_rescaled if time<=`H', fcolor(purple%15) lcolor(gs13) lw(none)) ///
				   (scatter result1 date_rescaled if time<=`H', c(l) clp(l) ms(i) clc(black) mc(black) clw(medthick)), ///
				   xline(`irf_line', lcolor(black) lwidth(medium)) ///
				   text(`text_y' `text_x' "← IRF Begins", place(e) size(small)) ///
				   xlabel(`xlabs', format(%tmCY)) xtitle("") yscale(range(. `=`y_max' + `y_range'*0.08')) ///
				   legend(order(1 "Decomposition" 3 "IRF") position(6) ring(0) cols(1)) ///
				   title("Decomposition + IRF: `y' response to `x'")

```

To add another line, simply do something like 

```
(line makes7 date if !missing(makes7), lcolor(blue) lwidth(medthick)) ///
```



