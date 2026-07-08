from objective import *
from BreguetRangeComp import *
from aero import *
from total_mass_comp import *
from propAndCost import *
from weight import *
import openmdao.api as om
import unittest
from openmdao.utils.assert_utils import assert_check_partials

class TestAero(unittest.TestCase):

    #inherits mnethods like:
    #self.assertEqual() et cetera

    #I guess this framework uses this name convention
    #runs before every func that starts with test__
    def setUp(self):
        self.prob = om.Problem()
        
        #dummy model to test aero
        self.prob.model.add_subsystem('Aero',AeroDiscipline(),promotes=['*'])
        #promotes makes sure evrything is accesible at self level
        #runs after every, individual, function that starts with test__
        #i guess its used to reset state of object being tested if needed

        self.prob.setup(force_alloc_complex=True)
        
        self.prob.set_val('m_total', 73229.6)
        self.prob.set_val('g', 9.80665)
        self.prob.set_val('rho', 0.38)
        self.prob.set_val('V', 230.0)
        self.prob.set_val('S', 102.0)
        self.prob.set_val('AR', 9.0)

        self.prob.set_val('C_D0_base', 0.02)
        self.prob.set_val('ks_base', 0.0005)
        self.prob.set_val('S_0', 100.0)
        self.prob.set_val('e_base', 0.8)

        self.prob.set_val('delta_CD0', 1.0)
        self.prob.set_val('delta_ks', 1.0)
        self.prob.set_val('delta_e', 1.0)

        self.prob.run_model()

    def tearDown(self):
        pass #does nothing 

    def test_partials(self):
        partial_data = self.prob.check_partials(method='cs')
        assert_check_partials(partial_data, atol=1e-12, rtol=1e-12)

    def test_bahavior(self):
        m_total = self.prob.get_val('m_total')
        g = self.prob.get_val('g')
        rho = self.prob.get_val('rho')
        V = self.prob.get_val('V')
        S = self.prob.get_val('S')

        expected_CL = (m_total * g) / (0.5 * rho * V**2 * S)
        actual_CL = self.prob.get_val('CL')

        np.testing.assert_allclose(actual_CL, expected_CL, rtol=1e-12, atol=1e-12)

        AR = self.prob.get_val('AR')

        C_D0_base = self.prob.get_val('C_D0_base')
        S_0 = self.prob.get_val('S_0')
        ks_base = self.prob.get_val('ks_base')
        e_base = self.prob.get_val('e_base')

        delta_CD0 = self.prob.get_val('delta_CD0')
        delta_ks = self.prob.get_val('delta_ks')
        delta_e = self.prob.get_val('delta_e')

        CL = (m_total * g) / (0.5 * rho * V**2 * S)

        C_D0 = C_D0_base * delta_CD0 + ks_base * delta_ks * (S - S_0)

        expected_CD = C_D0 + (CL**2) / (np.pi * AR * e_base * delta_e)
        actual_CD = self.prob.get_val('CD')

        np.testing.assert_allclose(actual_CD, expected_CD, rtol=1e-12, atol=1e-12)

        LoD = self.prob.get_val('LD')

        expected_LoD = CL/expected_CD
        np.testing.assert_allclose(LoD, expected_LoD, rtol=1e-12, atol=1e-12)

