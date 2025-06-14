program define standshock
    syntax varname, [NOISe(real 0.01)]
    
    * Store input variable name and create output name
    local input_var `varlist'
    local output_var `input_var'_std
    
    * Verify the new variable name doesn't exist
    capture confirm new variable `output_var'
    if _rc {
        display as error "Variable `output_var' already exists"
        exit 110
    }
    
    * Step 1: Create variable with only non-zero elements
    tempvar nonzero_only
    qui gen `nonzero_only' = `input_var' if `input_var' != 0
    
    * Step 2: Compute mean and standard deviation of non-zero elements
    quietly sum `nonzero_only'
    local mean = r(mean)
    local sd = r(sd)
    
    * Check if there are non-zero values
    if r(N) == 0 {
        display as error "No non-zero values found in `input_var'"
        exit 198
    }
    
    * Step 3: Standardize non-zero elements
    tempvar standardized_nonzero
    qui gen `standardized_nonzero' = (`nonzero_only' - `mean') / `sd'
    
    * Step 4: Create new variable where zeros remain zeros
    tempvar standardized_temp
    qui gen `standardized_temp' = `input_var'
    qui replace `standardized_temp' = `standardized_nonzero' if `input_var' != 0
    
    * Step 5: Rescale the entire series to have std = 1
    quietly sum `standardized_temp'
    local final_sd = r(sd)
    
    if `final_sd' == 0 {
        display as error "Final standard deviation is zero"
        exit 198
    }
    
    tempvar scaled
    qui gen `scaled' = `standardized_temp' / `final_sd'
    
    * Step 6: Add noise to non-zero values (if noise parameter > 0)
    if `noise' > 0 {
        tempvar random_noise
        qui gen `random_noise' = rnormal(0, `noise')
        qui gen `output_var' = `scaled' + `random_noise' if !missing(`scaled') & abs(`scaled') > 0
        qui replace `output_var' = `scaled' if missing(`output_var')
    }
    else {
        qui gen `output_var' = `scaled'
    }
    

end
