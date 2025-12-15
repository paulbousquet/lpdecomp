*! version 1.1.0  
program lpdecomp, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], [H(integer 20) H1(integer 0) Lag(integer 0) NWLag(integer 0) vmat(string) contemp(varlist) irfscale(integer 1) adjstd(real 1) ztail(real .05) zerod(varlist) sdep(varlist) CUM NOADJ NOX MULT NODRAW DELAY FORCE DDY DECOMPgraph MAKEScol(integer 1)]
	
	foreach var of local sdep {
		quietly tab `var'
		if r(r) != 2 {
			di as error "States must be binary."
			exit 459
		}
	}
	
	local fresh result1 irc1 irc2 time makes1 makes_date date_rescaled
	foreach var of local fresh {
		capture drop `var'
	}
	
	// Capture tsset info before preserve
	qui tsset
	local timevar = r(timevar)
	local tsfmt : format `timevar'
	
	preserve
	
	local y : word 1 of `anything'
	local x 
	local trc 
	local H = `h'
	local h1 = `h1'
	local nlag = `h'
	
	local mesh: subinstr local anything "`y'" "", word
	local ivdum = strpos("`mesh'", "(") 
	local x : word 1 of `mesh'
	local contr : subinstr local mesh "`x'" "", word
	local varlist `contr' `x' `y'
	
	local Lag = `lag'
	local w `contr'
	local vmat = "`vmat'"
	
	if (`Lag' > 0) {
		if ("`nox'" != "") {
			local lagvlist `w' `y'
		}
		else {
			if ("`ddy'" != ""){
				tempvar l_`y' 
				qui gen l_`y' = L.`y'
				local lagvlist `w' `x' `l_`y''
			}
			else {
				local lagvlist `varlist'
			}
		}
		if ("`delay'" != ""){
			local w 
		}
		foreach var of local lagvlist {
			forv i = 1/`Lag' {
				tempvar t`var'_`i'
				quietly gen `t`var'_`i'' = L`i'.`var'
				local w `w' `t`var'_`i''	
			} 
		}
	}
	
	local w `w' `contemp'
	scalar delt = 0 
	local EV = 1 
	
	quietly drop if _n <= `Lag'
	quietly drop if missing(`y')
	local T = _N
	
	if ("`cum'" != "") {
		forvalues i=`h1'/`H' {
			quietly gen temp_`i' = .
			forvalues j=1/`=`T'-`i'' {
				quietly sum `y' in `j'/`=`j'+`i''
				quietly replace temp_`i' = r(sum) in `j'
			}
		}
	}
	else {
		forvalues i=`h1'/`H' {
			quietly gen temp_`i'= F`i'.`y'
		}
	}
	
	local shabang `x' `w' 
	foreach var of local shabang {
		quietly drop if missing(`var')
	}
	
	local zerod `zerod'
	foreach var of local zerod {
		qui capture drop if `var'==0
	}

	tempvar h_range
	quietly generate `h_range' = `h1' + _n - 1 if _n <= `H'+1-`h1'
	quietly bspline, xvar(`h_range') power(1) knots(`=`h1''(1)`=`H'') gen(bs)
	mkmat bs*, matrix(bs)
	matrix basis = bs[1..`=`H'+1-`h1'',1'...]
	drop bs*
	
	if ("`w'" != "") {
		tempname wmat
		mkmat `w', matrix(`wmat')
		matrix w = (J(_N, 1, 1) , `wmat')
	}
	else {
		matrix w = J(_N, 1, 1)
	}
	
	local satur 
	local nvars : word count `sdep'
	local ncombos = 1 
	local K : word count `w'
	
	if (("`w'" != "") & (`ivdum'==0)) {
		quietly reg `x' `w'
		scalar delt = e(rmse)
	}
	else {
		if ("`w'" == "") {
			qui sum `x'
			scalar delt = r(sd)
		}
	}
	
	if ("`noadj'" != ""){
		local delta = 1
	}
	else {
		local delta = `=delt'*`irfscale'
	}
	
	forvalues i=`h1'/`H' {
		qui reg `x' `w'
		qui predict temp_x_`i', res 
		qui reg temp_`i' `w'
		qui predict temp_y_`i', res 
	}
	
	// Store the time values used in estimation
	mkmat `timevar', matrix(tused)
	
	mkmat temp_y_*, matrix(yy)
	mkmat temp_x_*, matrix(x)
	drop temp_*

	local T = _N
	local HR = `H' + 1 - `h1'
	local TS = `T' * `HR'
	local XS = colsof(basis)*`EV'
	local back = `T'-`h1'
	
	if (`ivdum' > 0){
		mkmat `trc', matrix(xz)
	}
	else {
		matrix xz = (1,1)
	}
	
	mata: yy = st_matrix("yy")
	mata: x = st_matrix("x")
	mata: basis = st_matrix("basis")
	mata: twirl(`back',yy,x,basis,`h1',`H',`HR',`TS',`XS',`EV', `K')
	mata: X = st_matrix("X")
	mata: Xoz = st_matrix("Xoz")
	mata: msc = st_matrix("msc")
	mata: sescl = st_matrix("sescl")
	mata: Y = st_matrix("Y") 
	mata: P = st_matrix("P") 
	mata: IDX = st_matrix("IDX")
	mata: xz = st_matrix("xz")
	mata: sel = st_matrix("sel")
	mata: grab = st_matrix("grab")
	mata: ivtwirl(`ivdum',`back',xz,basis,X,sel,`TS',`XS',`HR',`EV')
	mata: ZX = st_matrix("ZX")
	mata: cvtwirl(`T',Y,yy,X,Xoz,grab,msc,sescl,P,basis,IDX,ZX,`h1',`H',`=TS',`XS',`delta',`EV',`nlag',`ztail',"`vmat'",`adjstd')
	
	restore 
	
	svmat double results, names(result)
	svmat double irc, names(irc)
	svmat double makes, names(makes)
	svmat double tused, names(tused)
	
	gen time = _n - 1
	
	// Create date variable for makes aligned to estimation sample
	qui gen makes_date = tused1 if !missing(makes`makescol')
	format makes_date `tsfmt'
	
	// Find date range for decomposition
	qui sum makes_date if !missing(makes`makescol')
	local last_main = r(max)
	local first_main = r(min)
	
	// Vertical line exactly at end of decomposition
	local irf_line = `last_main'
	
	// IRF starts exactly at irf_line (time=0 maps to last_main)
	qui gen date_rescaled = `last_main' + time * 10 if time <= `H'
	format date_rescaled `tsfmt'
	
	// Get y range for text placement
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
	
	if ("`nodraw'" == ""){	
		if ("`decompgraph'" != "") {
			// Combined decomposition + IRF plot
			twoway (line makes`makescol' makes_date if !missing(makes`makescol'), lcolor(blue) lwidth(medthick)) ///
				   (rarea irc1 irc2 date_rescaled if time<=`H', fcolor(purple%15) lcolor(gs13) lw(none)) ///
				   (scatter result1 date_rescaled if time<=`H', c(l) clp(l) ms(i) clc(black) mc(black) clw(medthick)), ///
				   xline(`irf_line', lcolor(black) lwidth(medium)) ///
				   text(`text_y' `text_x' "← IRF Begins", place(e) size(small)) ///
				   xlabel(`xlabs', format(%tmCY)) xtitle("") yscale(range(. `=`y_max' + `y_range'*0.08')) ///
				   legend(order(1 "Decomposition" 3 "IRF") position(6) ring(0) cols(1)) ///
				   title("Decomposition + IRF: `y' response to `x'")
		}
		else if ("`mult'" == "") {
			// Standard IRF plot
			tw (rarea irc1 irc2 time, fcolor(purple%15) lcolor(gs13) lw(none)) ///
			   (scatter result1 time, c(l) clp(l) ms(i) clc(black) mc(black) clw(medthick) legend(off)) ///
			   if time<=`H', title("IRF of `y' for shock to `x'") xtitle("horizon")
		}
		else {
			forvalues i=1/`EV' {
				local xx : word `i' of `endg'
				local shift = (`H'+1)*(`i'-1)
				tempvar sc_time 
				quietly gen `sc_time' = time - `shift'
				tw (rarea irc1 irc2 `sc_time', fcolor(purple%15) lcolor(gs13) lw(none)) ///
				   (scatter result1 `sc_time', c(l) clp(l) ms(i) clc(black) mc(black) clw(medthick) legend(off)) ///
				   if (`shift'<=time)&(time<=`=`i'*(`H'+1)-1'), title("IRF of `xx' for shock to `xx'") name("`xx'") xtitle(horizon)
			}
		}
	}
	
	// Clean up internal temp variables (keep makes_date, irf_xpos for user)
	cap drop tused1
end

mata: 

void function twirl( real scalar back,
                     real matrix yy,
		     real matrix x,
		     real matrix basis,
		     real scalar h1,
		     real scalar H, 
		     real scalar HR,
		     real scalar TS,
		     real scalar XS,
		     real scalar EV,
		     real scalar K)
{

	IDX = J(TS, 2, .)
	grab = J(back+h1,1,(1::HR))
	Y = J(TS, 1, .)
	Xb = J(TS, XS, .)
	II = I(HR)
	XSt = XS/EV
	base = J(1, HR, TS/HR )
	msc = TS/HR :- colsum(yy :== .)
	sescl = sqrt((msc:-2):/(msc:-(2+K)))
		
	for(t=1; t<=back; t++) {
		stt = (t-1)*HR + 1
		edd = t*HR
    
		IDX[|stt,1 \ edd,2|] = J(HR, 1, t), range(h1, H,1)
   
		Y[|stt,1 \ edd,1|] = yy[|t,1 \ t,HR|]'
		
		for(i=1; i<=EV; i++) {
			Xb[|stt,idb(i,XSt) \ edd,i*XSt|] = x[t,.]' :* basis
		}
	
	}

	sel = Y :!= .
	IDX = select(IDX, sel)
	Y = select(Y, sel)
	X = select(Xb, sel)
	grab = select(grab,sel)
	TS = length(Y)
	X_tau = J(1, XS, 1)
	st_matrix("X", X)
	st_matrix("Xoz",X_tau)
	st_matrix("msc", msc')
	st_matrix("sescl", sescl')
	st_matrix("Y", Y) 
	st_matrix("IDX", IDX)
	st_matrix("grab",grab)
	st_numscalar("TS",TS)
	st_matrix("sel",sel)
	P = J(cols(X), cols(X), 0)
    D = I(XSt)
	DD = D' * D
    P[1::XSt, 1::XSt] = DD
	for(i=2; i<=EV; i++) {
		stt = XSt*(i-1)+1
		edd = XSt*i 
		P[stt::edd,stt::edd] = DD 
	}
	st_matrix("P",P)
} 

void function cvtwirl( real scalar T,
                     real matrix Y,
					 real matrix yy,
		     real matrix X,
			 real matrix Xoz, 
			 real matrix grab, 
			 real matrix msc, 
			 real matrix sescl,
		     real matrix P,
		     real matrix basis,
		     real matrix IDX,
			 real matrix ZX, 
		     real scalar h1,
		     real scalar H, 
		     real scalar TS,
		     real scalar XS,
		     real scalar delta,
		     real scalar EV,
			 real scalar nlag,
			 real scalar ztail,
		     string vmat,
			 real scalar adjstd)
{
	
	linked = (H + 1)*EV
        results = J(linked, 1, 0)
        theta = J(cols(X), 1, 0)
        lambda_opt = 10^(-10)
	
	XX = quadcross(X, X)
	XY = quadcross(X, Y)
	
        A = XX + lambda_opt * rows(Y) * P
		bread = luinv(A)
        theta = bread * XY
		thetav = theta 
	beta = theta[1..XS, 1]
	results[|(h1+1),1 \ linked,1 |] = basis * beta * delta
        st_matrix("results", results)
		u = Y - ZX * thetav
        S = X :* (u * J(1, cols(X), 1))

	lagseq = 0::nlag
        V = quadcross(S, S)
	meat = V 
	if (vmat=="nw"){
		lagseq = 0::nlag
		weights = 1 :- lagseq :/ (nlag+1)
		for (i=1; i<=nlag; i++){
		Gammai = quadcross(S[(i+1)::rows(S),], S[1::(rows(S)-i),])
		GplusGprime = Gammai + Gammai'	
		V = V + weights[i+1] * GplusGprime
		}
		meat = V
	}
        V = bread * meat * bread
		V = basis * V[1::XS, 1::XS] * basis'
		 se = sescl :* sqrt(diagonal(V))
		 conf = J(rows(se), 2, .)
		 mu = basis * beta
        conf[,1] = mu :+ se * invnormal(ztail)
        conf[,2] = mu :+ se * invnormal(1-ztail)
		irc = J(rows(se)+1, 2, .)
        irc[(h1+1)::linked,] = conf * delta
	results[|(h1+1),1 \ linked,1 |] = mu * delta
	st_matrix("results", results)
        st_matrix("irc", irc)
	st_matrix("se", se)
	wvec = delta * Xoz * bread * X'
	plate = wvec'
	bleak = plate :* Y
	conv = yy 
	wvec = conv 
	for (j = 1; j <= XS; j++) {
    sp = msc[j,]
	sel = grab :== j
    mike = select(bleak, sel)
	jerry = select(plate, sel)
	yy[1::sp, j] = runningsum(mike, 0) 
	conv[1::sp, j] = mike 
	wvec[1::sp, j] = jerry 
	
}
	st_matrix("wvec", wvec)
	st_matrix("conv",conv)
    st_matrix("makes", yy)
}

void function ivtwirl( real scalar ivdum,
                     real scalar back,
		     real matrix xz,
		     real matrix basis,
			 real matrix X, 
			 real matrix sel,
			 real scalar TS, 
		     real scalar XS,
			 real scalar HR, 
		     real scalar EV) 
{
	if (ivdum > 0){
		Xb = J(TS, XS, 0)
		for(t=1; t<= back; t++){
			idx_beg = (t-1)*HR + 1
			idx_end = t*HR
			stack = basis*xz[t,1]
			for(i=2; i<=EV; i++){
				stack = stack, basis*xz[t,i]
			}
			Xb[|idx_beg,1 \ idx_end,XS|] = stack
		}
		Xb = select(Xb,sel)
		ZX = Xb, X[1..rows(X),(XS+1)..cols(X)]
		st_matrix("ZX", ZX)
	}
	else{
	   ZX = X	
	}
	st_matrix("ZX",ZX)
}


real scalar function idb(idx,blk){
	return ((idx-1)*blk+1)
}
end 
