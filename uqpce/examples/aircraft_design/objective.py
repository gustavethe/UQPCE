import openmdao.api as om
import numpy as np
from fixed import parameters

class DOC(om.ExplicitComponent):
    """
    Component for "DOCComp" box containing analytical derivatives
    """
    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']

        #Parameters
        self.add_input('Cf_base', units='USD/kg')
        self.add_input('C_time', units='USD/s')
        self.add_input('k_acq')
        self.add_input('C_eng_ref', units='USD')
        self.add_input('beta_base')

        #Global design variables
        self.add_input('SFC_tech', val=0., desc='SFC technology factor')
        self.add_input('V', units='m/s', desc='Cruise speed')

        #Local design variable
        self.add_input('R', units='m', desc='Breguet range', shape=(n,))
        
        #Solver state
        self.add_input('m_fuel', units='kg', desc='Fuel mass', shape=(n,)) 

        #Uncertainties
        self.add_input('delta_Cf', val=1.0, shape=(n,))
        self.add_input('delta_beta', val=1.0, shape=(n,))

        #Output
        self.add_output('DOC', units='USD', desc="Direct operating cost", shape=(n,))

    def setup_partials(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.declare_partials('DOC', ['V', 'SFC_tech', 'Cf_base', 'C_time', 'k_acq', 'C_eng_ref', 'beta_base'])
        self.declare_partials('DOC', ['R', 'm_fuel', 'delta_Cf', 'delta_beta'], rows=arange, cols=arange)

    def compute(self, inputs, outputs):
        """
        DOC = Cf_base * delta_Cf * m_fuel + C_time * (R / V) + k_acq * C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)
        """

        SFC_tech = inputs['SFC_tech']
        V = inputs['V']
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
        V = inputs['V']
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
        partials['DOC', 'V'] = -C_time * (R / V**2)
        partials['DOC', 'SFC_tech'] = k_acq * C_eng_ref * (beta_base * delta_beta)

        partials['DOC', 'Cf_base'] = delta_Cf * m_fuel
        partials['DOC', 'C_time'] = R / V
        partials['DOC', 'k_acq'] = C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)
        partials['DOC', 'C_eng_ref'] = k_acq * (1 + beta_base * delta_beta * SFC_tech)
        partials['DOC', 'beta_base'] = (k_acq * C_eng_ref) * (delta_beta * SFC_tech)

        partials['DOC', 'delta_Cf'] = Cf_base * m_fuel
        partials['DOC', 'delta_beta'] = (k_acq * C_eng_ref) * (beta_base * SFC_tech)

        # partials['Dpm', 'm_fuel'] = partials['DOC', 'm_fuel'] / (N_pax * R)
        # partials['Dpm', 'R'] = -(Cf_base * delta_Cf * m_fuel + k_acq * C_eng_ref * (1 + beta_base * delta_beta * SFC_tech)) / (N_pax * R**2)
        # partials['Dpm', 'V'] = partials['DOC', 'V'] / (N_pax * R)
        # partials['Dpm', 'SFC_tech'] = partials['DOC', 'SFC_tech'] / (N_pax * R)
        # partials['Dpm', 'Cf_base'] = partials['DOC', 'Cf_base'] / (N_pax * R)
        # partials['Dpm', 'C_time'] = partials['DOC', 'C_time'] / (N_pax * R)
        # partials['Dpm', 'k_acq'] = partials['DOC', 'k_acq'] / (N_pax * R)
        # partials['Dpm', 'C_eng_ref'] = partials['DOC', 'C_eng_ref'] / (N_pax * R)
        # partials['Dpm', 'beta_base'] = partials['DOC', 'beta_base'] / (N_pax * R)
        # partials['Dpm', 'N_pax'] = -(DOC / (N_pax**2 * R))

        # partials['Dpm', 'delta_Cf'] = partials['DOC', 'delta_Cf'] / (N_pax * R)
        # partials['Dpm', 'delta_beta'] = partials['DOC', 'delta_beta'] / (N_pax * R)

class Dpm(om.ExplicitComponent):
    """
    Component for objective of minimizing DOC/pax*km
    """
    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']

        #Parameters
        self.add_input('DOC', units='USD', shape=(n,))

        self.add_input('N_pax',val=189, desc="Number of passengers")

        #Local design variable
        self.add_input('R', units='km', desc='Breguet range', shape=(n,))

        #Output
        self.add_output('Dpm', desc="DOC/pax*km", shape=(n,))

    def setup_partials(self):
        n = self.options['vec_size']
        arange = np.arange(n)

        self.declare_partials('Dpm', ['N_pax'])
        self.declare_partials('Dpm', ['R', 'DOC'], rows=arange, cols=arange)

    def compute(self, inputs, outputs):
        """
        Dpm = DOC / (pax * km)
        """

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



class DesignMatch737Objective(om.ExplicitComponent):
    """
    Objective that penalizes distance from 737-800-like design variables.

    Uses normalized squared error so S, AR, and V are comparable.
    """

    def setup(self):
        self.add_input('S', val=parameters['S'], units='m**2')
        self.add_input('AR', val=parameters['AR'])
        self.add_input('V', val=parameters['V'], units='m/s')
        self.add_input('SFC_tech', val=parameters['SFC_tech'])

        self.add_output('J_737', val=0.0)

    def setup_partials(self):
        self.declare_partials('J_737', ['S', 'AR', 'V', 'SFC_tech'])

    def compute(self, inputs, outputs):
        S = inputs['S']
        AR = inputs['AR']
        V = inputs['V']
        SFC_tech = inputs['SFC_tech']

        S_ref = parameters['S']
        AR_ref = parameters['AR']
        V_ref = parameters['V']

        outputs['J_737'] = (
            ((AR - AR_ref) / AR_ref) ** 2
            + ((S - S_ref) / S_ref) ** 2
            + ((V - V_ref) / V_ref) ** 2
            + SFC_tech ** 2
        )

    def compute_partials(self, inputs, partials):
        S = inputs['S']
        AR = inputs['AR']
        V = inputs['V']
        SFC_tech = inputs['SFC_tech']

        S_ref = parameters['S']
        AR_ref = parameters['AR']
        V_ref = parameters['V']

        partials['J_737', 'AR'] = 2.0 * (AR - AR_ref) / (AR_ref ** 2)
        partials['J_737', 'S'] = 2.0 * (S - S_ref) / (S_ref ** 2)
        partials['J_737', 'V'] = 2.0 * (V - V_ref) / (V_ref ** 2)
        partials['J_737', 'SFC_tech'] = 2.0 * SFC_tech