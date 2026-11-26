"""Small analytical checks; run locally, not on push."""
import unittest
import numpy as np
from pyaehmop.transport import principal_layers, spectrum

class TransportTest(unittest.TestCase):
    def test_uniform_chain(self):
        R=np.array([[-1],[0],[1]])
        H=np.array([[[-1.]],[[0.]],[[-1.]]])
        S=np.array([[[0.]],[[1.]],[[0.]]])
        t,_=spectrum(R,H,S,[-1.,.3,3.],[[]],eta=1e-5)
        np.testing.assert_allclose(t[:2],1.,atol=1e-3)
        self.assertLess(t[2],1e-6)
    def test_long_range_folding(self):
        R=np.array([[-2],[0],[2]])
        H=np.array([[[-1.]],[[0.]],[[-1.]]])
        S=np.array([[[0.]],[[1.]],[[0.]]])
        h,s=principal_layers(R,H,S,[])
        self.assertEqual(h.shape,(3,2,2))
        np.testing.assert_allclose(h[2],-np.eye(2))
        np.testing.assert_allclose(s[1],np.eye(2))

if __name__=="__main__": unittest.main()
