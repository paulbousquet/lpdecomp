*! version 1.0.0  
program flpdecomp, eclass 
    version 15.0
    syntax anything(equalok) [if] [in], [H(integer 20) H1(integer 0) Lag(integer 0) NWLag(integer 0) vmat(string) contemp(varlist) irfscale(integer 1) adjstd(real 1) ztail(real .05) zerod(varlist) sdep(varlist) CUM NOADJ NOX NODRAW DELAY FORCE]
	
	foreach var of local sdep {
		quietly tab `var'
		if r(r) != 2 {
			di as error "States must be binary. If you have multiple, discrete levels, this can also be coded as mulitple binary variables."
        exit 459
		}
	}
	
	local fresh result1 irc1 irc2 time 
	
	foreach var of local fresh {
    capture drop `var'
	}
	preserve
    // Parse input
    local y : word 1 of `anything'

    // Initialize x 
    local x 
	local trc 
    local H = `h'
    local h1 = `h1'
	local nlag = `h'
	
    // Remove the dependent variable from the variable list
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
			local lagvlist `varlist'
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
	
    	local EV  = 1 
    
    
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
	

    // Create the range of values for h
    tempvar h_range
    quietly generate `h_range' = `h1' + _n - 1 if _n <= `H'+1-`h1'

    // Generate basis functions
    quietly bspline, xvar(`h_range') power(1) knots(`=`h1''(1)`=`H'') gen(bs)
    
    mkmat bs*, matrix(bs)

    // Create a matrix from the basis variables
    matrix basis = bs[1..`=`H'+1-`h1'',1'...]
	    
    drop bs*
    
    // Create covariates matrix if w is specified
    
        if ("`w'" != "") {
        tempname wmat
        mkmat `w', matrix(`wmat')
        matrix w = (J(_N, 1, 1) , `wmat')
    }
    else {
        matrix w = J(_N, 1, 1)
    }
    
    // 1 shock std dev
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
	local satur 
	
	local nvars : word count `sdep'
	local ncombos = 2^`nvars'

	quietly forvalues i = 0/`=`ncombos'-1' {
		local combo_num = `i' + 1
		gen combo`combo_num' = 1
    
		forvalues j = 1/`nvars' {
			local bit = 1-mod(int(`i'/2^(`j'-1)), 2)
			local var : word `j' of `sdep'
			replace combo`combo_num' = combo`combo_num' & (`var' == `bit')
		}
		if (`i' != 1){
			local satur `satur' combo`combo_num'
		}
	}
	
	
		local ww `w'
		local xx `x'
		local w 
		foreach wvar of local ww {
			forvalues i=1/`ncombos'{
				tempvar `wvar'_`i'
				gen ``wvar'_`i'' = `wvar' * combo`i'
				local w `w' ``wvar'_`i''
			}
		}
		foreach wvar of local xx {
			tempvar `wvar'_1
			gen ``wvar'_1' = `wvar' * combo1 
			local x ``wvar'_1'
			forvalues i=2/`ncombos'{
				tempvar `wvar'_`i'
				gen ``wvar'_`i'' = `wvar' * combo`i'
				local w `w' ``wvar'_`i''
			}
		}

if ("`satur'" != "" & "`force'" != "") {
            foreach var of local satur {
                qui sum `var'
                if r(mean) < .01 {
                    di as error "One of your states covers less than 1% of the sample. Double check that there are no trivial inclusions. If you are comfortable with this choice, re-run this command with the 'force' option added."
                    exit 459
                }
            }
        }
local done 0
while !`done' {
    local missing 0
    reg `x' `satur' `w' if !missing(temp_`h1')
    
    if e(F) == . local missing 1
    di `missing'
    
    if `missing' {
        local n : word count `w'
        
        display "/******************************************************************************/"
        display "/* This specification likely suffers from the curse of dimensionality (small T/K) */"
        display "/*                                                                              */"
        display "/* Variables will be dropped based on the eigenvalues of the design matrix       */"
        display "/******************************************************************************/"
        display "Number of controls at start: `=`n''"
        
        qui _rmcoll `w'
        local w = r(varlist)
        local w : display ustrregexra("`w'", "o\.[^ ]*", "")
        
        qui correlate `w'
        matrix C = r(C)
        matrix symeigen X L = C
        matrix L = L'
        
        mata: st_matrix("drop_ind", (st_matrix("L") :<= 0.001))
        
        mata: {
            w_vars   = tokens(st_local("w"))
            drop_vars = select(w_vars, st_matrix("drop_ind")')
            st_local("todrop", invtokens(drop_vars))
            keep_vars = select(w_vars, !st_matrix("drop_ind")')
            st_local("w", invtokens(keep_vars))
        }
        
        local newn : word count `w'
        local drp = `n' - `newn'
        
        if `drp' == 0 {
            display "/******************************************************************************/"
            display "/* A more aggressive approach is needed for variable selection                 */"
            display "/******************************************************************************/"
            
            local keep = `n' - 10
            qui lasso linear `x' (`satur') `w' if !missing(temp_`h1')
            
            di e(k_allvars)
            local lank = e(lambda_last)
            local lank_prev = `lank'
            local trim = 0 
            
            while !`trim' {
                if (e(k_nonzero_max) < `keep' & `lank' > 1e-15) {
                    display `=`keep' - e(k_nonzero_max)'
                    local lank = .1 * `lank'
                    
                    capture lasso linear `x' (`satur') `w' if !missing(temp_`h1'), ///
                        lambda(`lank') selection(none)
                    
                    if (rc == 0) {
                        local lank_prev = `lank'
                    }
                }
                else {
                    local trim 1 
                }
            }
            
            qui lasso linear `x' (`satur') `w' if !missing(temp_`h1'), lambda(`lank_prev')
            local w = e(allvars_sel)
            
            local newn : word count `w'
            local drp = `n' - `newn'
            
            if `drp' == 0 {
                di as error "Automatic variable selection failed. You must manually reduce the control set."
                exit 459
            }
            else {
                if (`lank_prev' != `lank') {
                    display "/******************************************************************************/"
                    display "/* WARNING: there were convergence issues. Number of variables dropped may be significant. */"
                    display "/******************************************************************************/"
                    display "`=`drp'' controls dropped"
                }
                else {
                    display "`=`drp'' controls dropped"
                }
            }
        }
        else {
            display "`=`drp'' controls dropped"
        }
    }
    else {
        local done 1
    }
}

		forvalues i=`h1'/`H' {
		 reg `x' `w' if !missing(temp_`i')
		qui predict temp_x_`i', residuals 
		qui reg temp_`i' `w'
		qui predict temp_y_`i',residuals 
		}
		mkmat temp_y_*, matrix(yy)
		mkmat temp_x_*, matrix(x)
		drop temp_*
	
    // Additional initializations
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
    mata: twirl(`back',yy,x,basis,`h1',`H',`HR',`TS',`XS',`EV')
	mata: X = st_matrix("X")
	mata: Xoz = st_matrix("X")
	mata: msc = st_matrix("X")
    mata: Y = st_matrix("Y") 
    mata: P = st_matrix("P") 
    mata: IDX = st_matrix("IDX")
	mata: sel = st_matrix("sel")
	mata: ZX = st_matrix("ZX")
	mata: grab = st_matrix("X")
    mata: cvtwirl(`T',Y,yy,X,Xoz,grab,msc,P,basis,IDX,ZX,`h1',`H',`=TS',`XS',`delta',`EV',`nlag',`ztail',"`vmat'",`adjstd')
	restore 
        // Prepare data for graphing
    svmat double results, names(result)

        
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
		     real scalar EV)
{

	IDX = J(TS, 2, .)
	grab = J(back+h1,1,(1::HR))
	Y = J(TS, 1, .)
	Xb = J(TS, XS, .)
	II = I(HR)
	XSt = XS/EV
	
		
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

	
	st_matrix("X", X)
	st_matrix("Y", Y) 
	st_matrix("IDX", IDX)
	st_numscalar("TS",TS)
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


}


real scalar function idb(idx,blk){
	return ((idx-1)*blk+1)
}
end 