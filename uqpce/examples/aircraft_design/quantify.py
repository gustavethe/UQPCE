import openmdao.api as om
import numpy as np
import matplotlib.pyplot as plt
#from disciplines.objective import *
from disciplines.BreguetRange import *
from disciplines.aero import *
from disciplines.total_mass_comp import *
from disciplines.propAndCost import *
from disciplines.weight import *
from disciplines.objective import *

from helpers import *
from sweepers import *



class CoupledDisciplines(om.Group):

    def initialize(self):
        self.options.declare('vec_size',default=1, types=int)

    def setup(self):
        n = self.options['vec_size']

        ###Total Mass Component####################################
        self.add_subsystem(
            'Mass',TotalMassComp(vec_size=n),
            promotes_inputs=['m_empty',
                             'm_fuel'],
            promotes_outputs=['m_total']
                           )
        #^######################################################^#

        ###Breguet Range Component################################
        self.add_subsystem(
            'Range',BreguetRangeComp(vec_size=n),
            promotes_inputs=['V_cruise',
                            'm_total','LD',
                            'SFC',
                            'm_fuel'],
            promotes_outputs=['R']
                           )
        #^######################################################^#

        ###Structural Weight Component############################
        self.add_subsystem(
            'Weight',Weights_Struct(vec_size=n),
            promotes_inputs=['S','AR','V_cruise',
            'delta_kw','delta_fsys','delta_p',
            'm_total','m_engine'],
            promotes_outputs=['m_wing','m_empty']
                           )
        #^######################################################^#

        ###Aerodynamics Component#################################
        self.add_subsystem(
            'Aero',AeroComp(vec_size=n),
            promotes_inputs=['S','AR','V_cruise',
            'delta_CD0','delta_ks','delta_e',
            'm_total'], 
            promotes_outputs=['CL','CD','LD','WL']
                           )
        #^######################################################^#

        ###Range Residual#########################################
        initial_guess = np.ones(n)*16000 #kg
        Balance = om.BalanceComp()
        
        Balance.add_balance(
            name='m_fuel',val=initial_guess,
            units='kg',lower=1000.0,upper=50000.0,
            lhs_name='R',rhs_name='R_target',
            rhs_val=parameters['R_target'],
            eq_units='m',ref=16000.0,res_ref=1.0e6,
            )
        
        self.add_subsystem('Balance', Balance,
                           promotes_inputs=['R'],
                           promotes_outputs=['m_fuel'])
        #^######################################################^#
        
        ###Residual Solver Options################################
        newton = self.nonlinear_solver = om.NewtonSolver(solve_subsystems=True)
        self.nonlinear_solver.options['iprint'] = 2
        self.nonlinear_solver.options['maxiter'] = 500
        self.nonlinear_solver.options['atol'] = 1e-7
        self.nonlinear_solver.options['rtol'] = 1e-7

        line_search = newton.linesearch = om.ArmijoGoldsteinLS(
                                    bound_enforcement='vector',
                                        )
        line_search.options['maxiter'] = 20
        line_search.options['print_bound_enforce'] = True
        self.linear_solver = om.DirectSolver()
        #^######################################################^#

class ExampleMDA(om.Group):

    def initialize(self):
        self.options.declare('vec_size',default=1, types=int)
    
    def setup(self):
        n = self.options['vec_size']

        ###Propulsion Components##################################
        self.add_subsystem(
            'Prop', Propulsion(vec_size=n),
            promotes_inputs=['delta_eta','delta_kv',
                             'V_cruise','SFC_tech'],
            promotes_outputs=['SFC']
                          )
        #^######################################################^#

        ###Engine Weight Component################################
        self.add_subsystem(
            'Engine', EngineWeight(vec_size=n), 
            promotes_inputs=['delta_alpha','SFC_tech'],
            promotes_outputs=['m_engine']
                           )
        #^######################################################^#
        
        ###Coupled Component Group################################
        self.add_subsystem(
            'Coupled', CoupledDisciplines(vec_size=n), 
            promotes_inputs=['delta_kw','delta_fsys','delta_p',
                            'delta_CD0','delta_ks','delta_e',
                            'S','AR','V_cruise',
                            'SFC','m_engine'],
            promotes_outputs=['R',
                              'm_fuel','m_total',
                              'm_empty','m_wing',
                              'CL','CD','LD','WL']
                           )
        #^######################################################^#

class DOC(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']

        #proposed design variables
        self.add_input('SFC_tech', units=None)
        self.add_input('V_cruise', units='m/s')
        
        #model variable (output from other component)
        self.add_input('R', units='m', shape=(n,))
        self.add_input('m_fuel', units='kg',shape=(n,)) 

        #uncertain parameters
        self.add_input('delta_Cf',val=np.ones(n),units=None,shape=(n,))
        self.add_input('delta_beta',val=np.ones(n), units=None,shape=(n,))

        #tuning parameters
        self.add_input('Cf_base', units='USD/kg')
        self.add_input('beta_base', units=None)
        
        #constant parameters
        self.add_input('C_time', val=parameters['C_time'], units='USD/s')
        self.add_input('k_acq', val=parameters['k_acq'], units=None)
        self.add_input('C_eng_ref', val=parameters['C_eng_ref'], units='USD')

        #outputs
        self.add_output('DOC', units='USD',shape=(n,))
       
    def setup_partials(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.declare_partials('DOC', ['V_cruise', 'SFC_tech', 'Cf_base', 'C_time', 'k_acq', 'C_eng_ref', 'beta_base'])
        self.declare_partials('DOC', ['R', 'm_fuel', 'delta_Cf', 'delta_beta'], rows=arange, cols=arange)

    def compute(self, inputs, outputs):
 
        SFC_tech = inputs['SFC_tech']
        V = inputs['V_cruise']
        Cf_base = inputs['Cf_base']
        m_fuel = inputs['m_fuel']
        C_time = inputs['C_time']
        R = inputs['R']
        k_acq = inputs['k_acq']
        C_eng_ref = inputs['C_eng_ref']
        beta_base = inputs['beta_base']
        delta_beta = inputs['delta_beta']
        delta_Cf = inputs['delta_Cf']

        outputs['DOC'] = DOC = Cf_base * delta_Cf * m_fuel + C_time * (R/V) + k_acq * C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)
        
    def compute_partials(self, inputs, partials):
        SFC_tech = inputs['SFC_tech']
        V = inputs['V_cruise']
        Cf_base = inputs['Cf_base']
        m_fuel = inputs['m_fuel']
        C_time = inputs['C_time']
        R = inputs['R']
        k_acq = inputs['k_acq']
        C_eng_ref = inputs['C_eng_ref']
        beta_base = inputs['beta_base']
        delta_Cf = inputs['delta_Cf']
        delta_beta = inputs['delta_beta']

        # DOC = Cf_base * delta_Cf * m_fuel + C_time * (R/V) + k_acq * C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)

        partials['DOC', 'm_fuel'] = Cf_base * delta_Cf
        partials['DOC', 'R'] = C_time / V
        partials['DOC', 'V_cruise'] = -C_time * (R / V**2)
        partials['DOC', 'SFC_tech'] = k_acq * C_eng_ref * (beta_base * delta_beta)

        partials['DOC', 'Cf_base'] = delta_Cf * m_fuel
        partials['DOC', 'C_time'] = R / V
        partials['DOC', 'k_acq'] = C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)
        partials['DOC', 'C_eng_ref'] = k_acq * (1 + beta_base * delta_beta * SFC_tech)
        partials['DOC', 'beta_base'] = (k_acq * C_eng_ref) * (delta_beta * SFC_tech)

        partials['DOC', 'delta_Cf'] = Cf_base * m_fuel
        partials['DOC', 'delta_beta'] = (k_acq * C_eng_ref) * (beta_base * SFC_tech)

class Dpm(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']

        #proposed design variables
        #n/a

        #model variable (output from other component)
        self.add_input('DOC', units='USD', shape=(n,))
        self.add_input('R', units='km',shape=(n,))

        #uncertain parameters
        #n/a

        #tuning parameters
        #n/a

        #constant parameters
        self.add_input('N_pax', val=parameters['N_pax'])

        #outputs
        self.add_output('Dpm', shape=(n,))

    def setup_partials(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.declare_partials('Dpm', ['N_pax'])
        self.declare_partials('Dpm', ['R', 'DOC'], rows=arange, cols=arange)

    def compute(self, inputs, outputs):

        N_pax = inputs['N_pax']
        DOC = inputs['DOC']
        R = inputs['R']

        outputs['Dpm'] = DOC / (N_pax * R)
    
    def compute_partials(self, inputs, partials):
        N_pax = inputs['N_pax']
        DOC = inputs['DOC']
        R = inputs['R']

        partials['Dpm', 'R'] = -(DOC / (N_pax * R**2))
        partials['Dpm', 'N_pax'] = -(DOC / (N_pax**2 * R))
        partials['Dpm', 'DOC'] = 1 / (N_pax * R)
        
class CL_constraint(om.ExplicitComponent):
    
    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.add_input('CL', shape=(n,))

        self.add_input('CL_target', val=0.53)

        self.add_output('CL_constraint', shape=(n,))

        # should be identity matrix
        self.declare_partials('CL_constraint', 'CL', rows=arange, cols=arange)

    def compute(self, inputs, outputs):

        CL = inputs['CL']
        CL_target = inputs['CL_target']

        outputs['CL_constraint'] = CL_target - CL

    def compute_partials(self, inputs, partials):

        partials['CL_constraint', 'CL'] = -1

class WingLoad_constraint(om.ExplicitComponent):
    
    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.add_input('WL', shape=(n,))

        self.add_input('WL_target',val=5905.0)

        self.add_output('WL_constraint', shape=(n,))

        # should be identity matrix
        self.declare_partials('WL_constraint', 'WL', rows=arange, cols=arange)

    def compute(self, inputs, outputs):

        WL = inputs['WL']
        WL_target = inputs['WL_target']

        outputs['WL_constraint'] = WL - WL_target

    def compute_partials(self, inputs, partials):

        partials['WL_constraint', 'WL'] = 1


from uqpce.mdao.uqpcegroup import UQPCEGroup
from uqpce.mdao import interface
import os
from fixed import optimal

def uqpce_main_script():
    #---------------------------------------------------------------------------
    #                               Input Files
    #---------------------------------------------------------------------------

    script_dir = os.path.dirname(os.path.abspath(__file__))
    relative_yaml = 'input.yaml'
    relative_matrix = 'run_matrix_generated.dat'
    input_file = os.path.join(script_dir, relative_yaml)
    matrix_file  = os.path.join(script_dir, relative_matrix)

    #---------------------------------------------------------------------------
    #                   Setting up for UQPCE and design under uncertainty
    #---------------------------------------------------------------------------

    (
        var_basis, norm_sq, resampled_var_basis, 
        aleatory_cnt, epistemic_cnt, resp_cnt, order, variables, 
        sig, run_matrix
    ) = interface.initialize(input_file, matrix_file)
    
    prob = om.Problem()
    
    #---------------------------------------------------------------------------
    #                   Add Subsystems to Problem
    #---------------------------------------------------------------------------
    
    prob.model.add_subsystem(
        'MDA', 
        ExampleMDA(vec_size=resp_cnt), 
        promotes_inputs=(['V_cruise', 'S', 'AR', 'SFC_tech',
                          'delta_eta', 'delta_kv','delta_alpha',
                          'delta_CD0','delta_ks','delta_e',
                          'delta_fsys','delta_kw','delta_p']), 
        promotes_outputs=['m_fuel','m_empty','m_engine',
                          'm_total','CL','CD','WL','SFC','R']
    )


    prob.model.add_subsystem(
        'WingLoad_constraint', 
        WingLoad_constraint(vec_size=resp_cnt), 
        promotes_inputs=['WL'], 
        promotes_outputs=['WL_constraint']
    )

    prob.model.add_subsystem(
        'LiftCoeff_constraint', 
        CL_constraint(vec_size=resp_cnt), 
        promotes_inputs=['CL'], 
        promotes_outputs=['CL_constraint']
    )


    prob.model.add_subsystem(
        'DOC_objective', 
        DOC(vec_size=resp_cnt), 
        promotes_inputs=(['V_cruise','SFC_tech',
                          'delta_beta','delta_Cf','R','m_fuel']), 
        promotes_outputs=['DOC']
    )

    prob.model.add_subsystem(
        'DPM_objective', 
        Dpm(vec_size=resp_cnt), 
        promotes_inputs=['DOC','R'], 
        promotes_outputs=['Dpm']
    )



    #---------------------------------------------------------------------------
    #                   Add UQPCE Group to Problem
    #---------------------------------------------------------------------------

    probailistic_DOC_list = ['DOC:resampled_responses','DOC:ci_lower',
                             'DOC:ci_upper','DOC:mean','DOC:mean_plus_var']
    
    probailistic_m_fuel_list = ['m_fuel:resampled_responses','m_fuel:ci_lower',
                                'm_fuel:ci_upper','m_fuel:mean','m_fuel:mean_plus_var',]
    
    probailistic_m_empty_list = ['m_empty:resampled_responses','m_empty:ci_lower',
                                 'm_empty:ci_upper', 'm_empty:mean','m_empty:mean_plus_var',]
    
    probailistic_m_engine_list = ['m_engine:resampled_responses','m_engine:ci_lower',
                                  'm_engine:ci_upper','m_engine:mean','m_engine:mean_plus_var',]
    
    probailistic_m_total_list = ['m_total:resampled_responses','m_total:ci_lower',
                                 'm_total:ci_upper','m_total:mean','m_total:mean_plus_var',]
    
    probailistic_CL_list = ['CL:resampled_responses','CL:ci_lower',
                            'CL:ci_upper','CL:mean','CL:mean_plus_var']

    probailistic_CD_list = ['CD:resampled_responses','CD:ci_lower',
                            'CD:ci_upper','CD:mean','CD:mean_plus_var']
    
    probailistic_SFC_list = ['SFC:resampled_responses','SFC:ci_lower',
                             'SFC:ci_upper','SFC:mean','SFC:mean_plus_var',]
    
    probailistic_CL_constr_list = ['CL_constraint:resampled_responses',
                                   'CL_constraint:ci_lower',
                                   'CL_constraint:ci_upper',
                                   'CL_constraint:mean',
                                   'CL_constraint:mean_plus_var']

    probailistic_output_list = (probailistic_DOC_list +
                                probailistic_m_fuel_list +
                                probailistic_m_empty_list +
                                probailistic_m_engine_list +
                                probailistic_m_total_list +
                                probailistic_CL_list +
                                probailistic_CD_list +
                                probailistic_SFC_list +
                                probailistic_CL_constr_list)

    prob.model.add_subsystem(
        'UQPCE',
        UQPCEGroup(
            significance=sig,
            var_basis=var_basis,
            norm_sq=norm_sq,
            resampled_var_basis=resampled_var_basis,
            tail='both',
            epistemic_cnt=epistemic_cnt,
            aleatory_cnt=aleatory_cnt,
            uncert_list=['DOC', 'm_fuel','m_empty','m_engine','m_total','CL','CD','SFC','CL_constraint'],
            tanh_omega=1e-3,
            sample_ref0=[ 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0],
            sample_ref=[ 5.0e4, 1000, 1000, 1000, 1000,0.1,0.1,0.1,0.1],
        ),
        promotes_inputs=[ 'DOC', 'm_fuel','m_empty','m_engine','m_total','CL','CD','SFC','CL_constraint'],
        promotes_outputs= probailistic_output_list
    )


 


    #Assign objective function and constraint in UQPCE formatting
    
    #CL_con = 'CL_constraint:ci_lower'
    

    prob.model.set_input_defaults('S', val=optimal['S'], units='m**2')
    prob.model.set_input_defaults('AR', val=optimal['AR'])
    prob.model.set_input_defaults('V_cruise', val=optimal['V'], units='m/s')
    prob.model.set_input_defaults('SFC_tech', val=optimal['SFC_tech'])

    prob.driver = om.pyOptSparseDriver(optimizer='SLSQP')

    prob.model.add_design_var('S',lower=100.0,upper=180.0,ref=124.6,)

    prob.model.add_design_var('AR',lower=7.0,upper=50.0,ref=9.45)

    prob.model.add_design_var('V_cruise',lower=200.0,upper=260.0,ref=230.0)

    prob.model.add_design_var('SFC_tech',lower=-1.0,upper=1.0,ref=1.0)

    prob.model.add_objective('DOC:mean',ref=2.0e4)

    
    prob.model.add_constraint('CL_constraint:ci_upper',lower=0.0, ref0=1, ref=2)

    #prob.model.add_constraint('CL_constraint:ci_upper',upper=0.5, ref0=1, ref=2)


    #prob.model.add_constraint(
    #    'CL:mean',
    #    upper=0.53
    #)


    prob.setup()

    initialize(prob)
    interface.set_vals(
    prob,
    variables,
    run_matrix,
)

    prob.run_model()

    DOC_dist = prob.get_val('DOC:resampled_responses').copy().ravel()
    DOC_ci_lower = prob.get_val('DOC:ci_lower').copy().item()
    DOC_ci_upper = prob.get_val('DOC:ci_upper').copy().item()
    DOC_mu = prob.get_val('DOC:mean').copy().item()
    DOC_var_plus_mu = prob.get_val('DOC:mean_plus_var').copy().item()
    DOC_var = DOC_var_plus_mu - DOC_mu

    prob.check_totals(of=['DOC:mean','CL:mean',
        'CL:ci_lower',
        'CL:ci_upper',],wrt=['S', 'AR', 'SFC_tech','V_cruise'],
                      compact_print=True, method='fd')
    #initialize(prob)
    prob.run_driver()

    DOC_opt_dist = prob.get_val('DOC:resampled_responses').ravel()
    DOC_opt_ci_lower = prob.get_val('DOC:ci_lower').item()
    DOC_opt_ci_upper = prob.get_val('DOC:ci_upper').item()
    DOC_opt_mu = prob.get_val('DOC:mean').item()
    DOC_opt_var_plus_mu = prob.get_val('DOC:mean_plus_var').item()
    DOC_opt_var = DOC_opt_var_plus_mu - DOC_opt_mu

    
    fig, ax = plt.subplots()

    #fig.suptitle(r"Direct Operating Cost PDFs")

    ax.hist(DOC_dist,bins=100,density=True,color='purple',alpha=0.5)
    ax.axvline(DOC_ci_lower, color='red', linewidth=2,linestyle=':', label=rf"CI lower $\approx$ {DOC_ci_lower:.4f}")
    ax.axvline(DOC_ci_upper, color='red', linewidth=2,linestyle=':', label=rf"CI upper $\approx$ {DOC_ci_upper:.4f}")
    #ax.set_xlabel(r"$\mathrm{DOC}$ [USD]",labelpad=15,fontsize=18)
    #ax.set_ylabel(r"Probability Density",labelpad=10,fontsize=18)
    #ax.set_title(rf"Estimated DOC Distribution: $\mu = {DOC_mu:.4f}, \ \ \sigma^2 = {DOC_var:.4e}$",fontsize=24)
    
    ax.hist(DOC_opt_dist,bins=100,density=True,color='green',alpha=0.5)
    ax.axvline(DOC_opt_ci_lower, color='blue', linewidth=2,linestyle=':', label=rf"CI lower $\approx$ {DOC_ci_lower:.4f}")
    ax.axvline(DOC_opt_ci_upper, color='blue', linewidth=2,linestyle=':', label=rf"CI upper $\approx$ {DOC_ci_upper:.4f}")
    #ax.set_xlabel(r"$\mathrm{DOC}$ [USD]",labelpad=15,fontsize=18)
    #ax.set_ylabel(r"Probability Density",labelpad=10,fontsize=18)
    #ax.set_title(rf"Estimated DOC Distribution: $\mu = {DOC_mu:.4f}, \ \ \sigma^2 = {DOC_var:.4e}$",fontsize=24)
    ax.legend()


    plt.show()


    
   
 

def main():
    uqpce_main_script()

if __name__ == "__main__":
    main()