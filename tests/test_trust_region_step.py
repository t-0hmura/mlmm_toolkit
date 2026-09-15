"""Independent optimality and coordinate-space contracts for local fixes."""
import numpy as np
import pytest
import torch
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

@pytest.mark.parametrize('backend',['numpy','torch'])
@pytest.mark.parametrize('lam,g',[
    ([-1.,-1.,2.],[0.,0.,0.]),
    ([-1.,-1.,2.],[0.,0.,.03]),
    ([-1.,-1.,2.],[0.,1e-12,.03]),
    ([-1.,-1.+1e-10,2.],[0.,1e-7,0.]),
    ([0.,0.,2.],[0.,0.,.03]),
    ([1.,2.,3.],[.01,.02,.03]),
    ([1.,2.,3.],[1.,2.,3.]),
])
@pytest.mark.parametrize('transform',[False,True])
def test_trust_subproblem_kkt(lam,g,backend,transform):
    opt=RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius=.1; opt.log=lambda *_:None
    lam=np.array(lam);g=np.array(g)
    v=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
    physical_g=v@g
    args=[lam,v,physical_g]
    if backend=='torch': args=[torch.tensor(x,dtype=torch.float64) for x in args]
    with np.errstate(all='raise'):
        result=opt.get_newton_step_on_trust(*args,transform=transform)
    s=v.T@result if transform else result
    assert np.isfinite(s).all()
    assert np.linalg.norm(s)<=opt.trust_radius*(1+1e-12)
    if np.linalg.norm(s)<opt.trust_radius*(1-1e-10): shift=0.
    else: shift=-np.dot(s,lam*s+g)/np.dot(s,s)
    assert shift>=-1e-12
    assert np.min(lam+shift)>=-1e-10
    np.testing.assert_allclose((lam+shift)*s+g,0.,atol=1e-10)
    assert np.dot(s,g)+.5*np.dot(s,lam*s)<=1e-12

def test_unordered_lowest_eigenspace():
    opt=RFOptimizer.__new__(RFOptimizer);opt.trust_radius=.1;opt.log=lambda *_:None
    s=opt.get_newton_step_on_trust(np.array([2.,-1.,-1.]),np.eye(3),np.zeros(3))
    assert s[0]==0.
    assert np.linalg.norm(s)==pytest.approx(.1)

@pytest.mark.parametrize('seed',range(32))
def test_scalar_negative_hessian_bracket(seed):
    opt=RFOptimizer.__new__(RFOptimizer);opt.trust_radius=.1;opt.log=lambda *_:None
    g=np.random.default_rng(seed).normal(size=3)
    s=opt.get_newton_step_on_trust(-np.ones(3),np.eye(3),g)
    np.testing.assert_allclose(s,-.1*g/np.linalg.norm(g),rtol=1e-12,atol=1e-15)
