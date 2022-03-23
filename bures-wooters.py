import qutip
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pickle


def prod_basis(b1, b2):
    return [qutip.tensor(b,s) for b in b1 for s in b2]

def scalar_prod(op1,op2,rho0=None):
    if op1.dims[0][0]!=op1.dims[0][0]:
        return None
    if rho0 is None:
        rho0 = qutip.qeye(op1.dims[0])/op1.dims[0][0]
    return ((op1.dag()*op2+op2.dag()*op1)*rho0).tr()
        

def base_orto(ops,rho0):
    dim = ops[0].dims[0][0]
    base = []
    # hacer gramm schmidt
    for op in ops:
        coeff = [scalar_prod(op2,op, rho0) for op2 in base]
        op_mod = op - sum([ c*op2 for c, op2 in zip(coeff, base)])
        op_mod = op_mod/np.sqrt(scalar_prod(op_mod,op_mod, rho0))
        base.append(op_mod)
    return base


def proj_op(K,base,rho0):
    return sum([   scalar_prod(b, K,rho0) * b for b in base])

def logM(rho):
    vals, vecs = rho.eigenstates()
    return sum([np.log(val)*vec*vec.dag() for val,vec in zip(vals, vecs) if val>0])

def sqrtM(rho):
    vals, vecs = rho.eigenstates()
    return sum([ (abs(val)**.5)*vec*vec.dag() for val,vec in zip(vals, vecs)])

def rel_entropy(rho, sigma):
    val = (rho*(logM(rho)-logM(sigma))).tr()
    if abs(val.imag)>1.e-6:
        print("rho or sigma not positive")
        print(rho.eigenstates())
        print(sigma.eigenstates())
    return val.real


def bures(rho, sigma):
    val = abs((sqrtM(rho)*sqrtM(sigma)).tr())
    val = max(min(val,1.),-1.)
    return np.arccos(val)/np.pi
        
def maxent_rho(rho, basis):   
    def test(x, rho, basis):
        k = sum([-u*b for u,b in zip(x, basis)])        
        sigma = (.5*(k+k.dag())).expm()
        sigma = sigma/sigma.tr()
        return rel_entropy(rho, sigma)    
    res = opt.minimize(test,np.zeros(len(basis)),args=(rho,basis))
    k = sum([-u*b for u,b in zip(res.x, basis)])        
    sigma = (.5*(k+k.dag())).expm()
    sigma = sigma/sigma.tr()
    return sigma
    
    
def error_maxent_state(rho, basis, distance=bures):
    try:
        sigma = maxent_rho(rho, basis)
        return distance(rho,sigma)
    except:
        print("fail")
        return None
    
    
def error_proj_state(rho, rho0, basis, distance=bures):
    try:
        basis = base_orto(basis, rho0)
        sigma = proj_op(logM(rho), basis, rho0).expm()
        return distance(rho, sigma)
    except:
        print("fail")
        return None
class Result(object):
    def __init__(self, ts=None, states=None):
        self.ts = ts
        self.states = states
        self.max_ent_app = None
        self.projrho0_app = None
        self.projrho_inst_app = None 
        
dim1=7
dim2=8

def simul(omega_bos=3, omega_s=3, temp=1, gaussian=False, deltat=10., tmax=500., distance=bures):
    basis_bos1 = [qutip.qeye(dim1), qutip.create(dim1),qutip.create(dim1).dag(),qutip.num(dim1)]
    H_bos1 = qutip.tensor(qutip.num(dim1), qutip.qeye(dim2))
    rho0 = qutip.tensor(-.5*(qutip.num(dim1)/temp), -.5*(qutip.num(dim2)/temp)).expm()
    rho0 = rho0/rho0.tr()
    # Base
    if gaussian:
        basis_bos2 = [qutip.qeye(dim2), qutip.create(dim2),qutip.destroy(dim2).dag(),qutip.num(dim2)]
        H_bos2 = qutip.tensor(qutip.qeye(dim1), qutip.num(dim2))
    else:
        basis_bos2 = [qutip.qeye(dim2), qutip.num(dim2)]
        H_bos2 = qutip.tensor(qutip.qeye(dim1), qutip.num(dim2))
        
    basis = base_orto(prod_basis(basis_bos1, basis_bos2), rho0)
    H0 = omega_bos * H_bos1 + omega_s * H_bos2
    Hi = qutip.tensor(qutip.create(dim1),qutip.destroy(dim2))+qutip.tensor(qutip.destroy(dim1),qutip.create(dim2))
    H=H0+0.05*Hi
    # Hamiltoniano    
    
    sampling = int(10*max(1,omega_bos, omega_s)*deltat)
    
    states = [rho0]
    rho = rho0    
    ts= [0]
    for i in range(int(tmax/deltat)):
        result = qutip.mesolve(H, states[-1], np.linspace(0,deltat, sampling),args={'omega_bos': omega_bos, 'omega_s': omega_s})
        states.append(result.states[-1])
        ts.append(deltat*i)
    #[0.05*qutip.tensor(qutip.destroy(dim1),qutip.qeye(dim2)),0.1*qutip.tensor(qutip.create(dim1),qutip.qeye(dim2))]
    result = Result(ts, states)
    result.times = ts
    result.states = states
    result.max_ent_app = np.array([error_maxent_state(rho, basis, distance) for rho in states])
    result.projrho0_app = np.array([error_proj_state(rho, rho0, basis,distance) for rho in states])
    result.projrho_inst_app = np.array([error_proj_state(rho, qutip.tensor(rho.ptrace([0]),rho.ptrace([1])), 
                                                         basis, distance) for rho in states])
    
    if gaussian:
        title = f" BW Dinámica cerrada gaussiana wb1={omega_bos} wb2={omega_s} dim1={dim1} dim2={dim2}"
    else:
        title = f" BW Dinámica cerrada no gaussiana wb1={omega_bos} wb2={omega_s} dim1={dim1} dim2={dim2}" 

    with open(title+".pkl","wb") as f:
        pickle.dump(result, f)
    return result, title

## Dinámica no Gaussiana, no resonante

result, title = simul(omega_bos=3., omega_s=np.sqrt(48), temp=1, gaussian=True, deltat=5., tmax=500., distance=bures)


plt.plot(result.times, result.max_ent_app, color="orange", label="max-ent")
plt.plot(result.times, result.projrho0_app, color="violet", label="proj rho0")
plt.plot(result.times, result.projrho_inst_app, color="crimson", label="proj rho(t)")
plt.xlabel("t[s]")
plt.ylabel("Arccos(F)")

plt.legend()
plt.title(title)
plt.savefig(title + f" dim1={dim1}dim2={dim2}.svg")
