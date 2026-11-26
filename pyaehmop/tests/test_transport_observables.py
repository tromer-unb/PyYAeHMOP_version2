import unittest
import numpy as np
from pyaehmop.transport_observables import G0, conductance, current_voltage

class LandauerTests(unittest.TestCase):
    def test_constant_transmission(self):
        e=np.linspace(-3,3,601); t=np.ones_like(e); v=np.array([-1.,0.,1.])
        for temp in (0,300):
            np.testing.assert_allclose(current_voltage(e,t,v,temp),G0*v,atol=1e-11)
        np.testing.assert_allclose(conductance(e,t),G0)
    def test_grid_coverage(self):
        with self.assertRaises(ValueError):
            current_voltage([-1,1],[1,1],[3],300)
