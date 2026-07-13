import openmdao.api as om
from disciplines import *

class CoupledDisciplines(om.Group):

    def setup(self):

        self.add_subsystem('Aero',AeroDiscipline())
