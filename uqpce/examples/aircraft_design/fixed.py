#sourced from http://www.b737.org.uk/techspecsdetailed.htm
parameters = {
    #~~~~miscelaneous parameters~~~~~~~~~~~~~~~~~~
    "R_target": 5.5e6,          
    "N_pax": 189,              
    "SFC_ref": 1.60e-4,        
    "V_ref": 231.5,            
    "S_naught": 124.58,        
    "CD0_base": 0.022,         
    "e_oswald_base": 0.80,     
    "m_fuse": 14518,
    "m_payload_design": 17955.0, 
    "m_payload_max": 20540,
    "m_fuel_max": 21000,
    "m_wing" : 6941, 
    "m_eng_ref" : 8602,
    "m_total" : 50000,
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    "wing_load" : 5905,
    "AR" : 9.45,
    "S" : 124.58, 
    "V" : 240.5,
    "SFC_tech" : 0,
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #$$$$$COST STUFF$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    "Cf_base": 0.74,          
    "C_time": 1700.0 / 3600.0, 
    "k_acq": 0.00142, 
    "C_eng_ref": 2.2e7,        
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #extra
    "b" : 34.32, 
    "g" : 9.81,
    "rho" : 0.38,
}
optimal = {
    "CL" : 0.53,
    "S" : 170.29555507,
    "AR" : 10.57091064,
    "V" :231.68801345,
    "SFC_tech" : 0.45445116
}
tuning = {  
    "p_base" : 7.5443750000000005,
    "eta_base" : 0.4393500000352975,
    "kv_base" : 701.05144999999999,
    "alpha_base" : 0.345000107456725,
    "beta_base" : 0.55,
    "ks_base" : 0.0002910700075464138,
    "fsys_base": 0.19357,
    "kw_base": 53.0,
}




