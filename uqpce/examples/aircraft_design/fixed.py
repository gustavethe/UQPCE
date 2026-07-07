#hi
#737-800 parameters
#sourced from http://www.b737.org.uk/techspecsdetailed.htm
parameters = {
    #~~~~miscelaneous parameters~~~~~~~~~~~~~~~~~~
    "R_target": 4.6e6,           # try 4000 to 5000 km
    "N_pax": 189,              # num passengers
    "SFC_ref": 1.60e-4,        # 1/s, 737-class cruise fuel consumption rate
    "V_ref": 231.5,            # m/s, 737-800-ish cruise speed
    "S_naught": 124.58,        #[square meters] 737-800 wing area
    "CD0_base": 0.022,         # clean cruise parasite drag ballpark
    "e_oswald_base": 0.80,     # reasonable
    #~~~~miscelaneous parameters~~~~~~~~~~~~~~~~~~
    
    #~~~~~tuning parameters~~~~~~~~~~~~~~~~~~~~~~~~ 
    #determined by inspection using 737-800 values
    """
    "fsys_base": 0.19357,
    "kw_base": 53.0,
    "p_base": 4.3,     #1 until further notice...
    "eta_base": 0.360,       # shitty guesstimate
    "kv_base": 31.5,        # shitty guesstimate
    "alpha_base": 0.20,    # shitty guesstimate
    "beta_base": 0.21,       # shitty guesstimate
    "ks_base": 7.8e-5,    # SHITIEST guesstimate
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """

    "fsys_base": -0.028819162235424547,
    "kw_base":  105.78857378198913,
    "p_base": 7.529033675833872,     #1 until further notice...
    "eta_base": 0.4641840834222323,       # shitty guesstimate
    "kv_base": 41.085658256916716,        # shitty guesstimate
    "alpha_base": -0.5258823112563431,    # shitty guesstimate
    "beta_base": 0.38018335486778737,       # shitty guesstimate
    "ks_base":  0.00023970372729668102,    # SHITIEST guesstimate


    #~~revised mass estimates for 737-800 [kg]~~~~~
    "m_fuse": 14518,
    "m_payload_design": 17955.0, 
    "m_payload_max": 20540,
    "m_fuel_max": 21000, #GUESS init-val or COMPARE to out-val
    "m_wing" : 6941, #COMPARE to out-val
    "m_eng_ref" : 8602, #from 2 CFM56-7 
    "m_total" : 50000,
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~design var initializers based on 737-800~~~~
    "AR" : 9.45,
    "S" : 124.58, #same as 737-800 area
    "V" : 231.5,
    "SFC_tech" : 0,
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #$$$$$COST STUFF$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    "Cf_base": 0.80,           # cost/kg of fuel
    "C_time": 2000.0 / 3600.0, # cost/s
    "k_acq": 0.00122, #highly senssitive but interval seems to shift up with range
    "C_eng_ref": 2.0e7,        # normalized/reference total engine acquisition cost
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    #extra
    "b" : 34.32, # [meters] 737-800 wing span, use to compare
    "g" : 9.81,
    "rho" : 0.38,
}

optimal = {
    "CL" : 0.53,
    "S" : 143.59364805,
    "AR" : 18.15096907,
    "V" : 240.88550103,
    "SFC_tech" : 0.60919866
}






