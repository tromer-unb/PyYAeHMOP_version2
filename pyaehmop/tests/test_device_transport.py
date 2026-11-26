import unittest
import numpy as np
from pyaehmop.basis import Atom,Shell
from pyaehmop.device_transport import Region,build_device,spectrum_device
from pyaehmop.transport import realspace,spectrum
from pyaehmop.input import Structure

def atom(x,l=0):
    return Atom("H" if l==0 else "C",np.array([0.,0.,x]),1,(Shell(l+1,l,-10.,(1.3,),(1.,)),))
def region(x,l=0,**kw): return Region([atom(x,l)],**kw)
class DeviceTest(unittest.TestCase):
    def setUp(self):
        self.v=np.empty((0,3));self.r=np.array([0,0,1.4]);self.e=np.linspace(-12,-8,15)
    def make(self,L,D,R,spinor=False):
        return build_device(L,D,R,self.v,self.r,self.r,cutoff=1.5,quadrature=(16,12,8),spinor=spinor)
    def calc(self,M,channels=False): return spectrum_device(M,self.e,[[]],channels=channels)[0]
    def test_pristine_matches_original(self):
        M=self.make(region(-1.4),region(0),region(1.4))
        a=self.calc(M)[:,0]
        st=Structure("chain",[atom(0)],1,self.r[None,:],(2,),[],0,[])
        R,H,S=realspace(st,1.5,(16,12,8))
        b=spectrum(R,H,S,self.e,[[]])[0]
        np.testing.assert_allclose(a,b,atol=1e-7)
        self.assertGreater(a.max(),.9)
    def test_spin_degeneracy(self):
        a=self.calc(self.make(region(-1.4),region(0),region(1.4)))[:,0]
        b=self.calc(self.make(region(-1.4),region(0),region(1.4),True),True)
        np.testing.assert_allclose(b[:,0],2*a,atol=1e-7)
        np.testing.assert_allclose(b[:,1],b[:,2],atol=1e-7)
    def test_exchange(self):
        b=self.calc(self.make(region(-1.4),region(0,exchange={"0":{"s":.6}}),region(1.4),True),True)
        self.assertGreater(np.max(abs(b[:,1]-b[:,2])),1e-4)
        np.testing.assert_allclose(b[:,0],b[:,1]+b[:,2],atol=1e-8)
    def test_soc_and_rectangular_interfaces(self):
        for exchange in ({},{"0":{"p":.3}}):
            b=self.calc(self.make(region(-1.4),region(0,1,soc={"C":{"p":.2}},exchange=exchange),region(1.4),True))
            self.assertTrue(np.isfinite(b).all())
    def test_invalid_layer(self):
        with self.assertRaises(ValueError):
            build_device(region(-1.4),region(0),region(1.4),self.v,self.r,self.r,cutoff=3)
