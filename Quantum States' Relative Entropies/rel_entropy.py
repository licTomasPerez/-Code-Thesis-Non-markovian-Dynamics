import qutip
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pickle
import datetime

dim1=10
dim2=10

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
        #print(rho.eigenstates())
        #print(sigma.eigenstates())
    return val.real


def bures_metric(rho, sigma):
    val = abs((sqrtM(rho)*sqrtM(sigma)).tr())
    val = max(min(val,1.),-1.)
    return np.arccos(val)/np.pi
        
bures = bures_metric

def maxent_rho(rho, basis):   
    def test(x, rho, basis):
        k = sum([-u*b for u,b in zip(x, basis)])        
        sigma = (.5*(k+k.dag())).expm()
        sigma = sigma/sigma.tr()
        return rel_entropy(rho, sigma)    
    # x0 = np.zeros(len(basis))
    K = logM(rho)
    rho0 = qutip.tensor(rho.ptrace([0]),rho.ptrace([1]))
    basis = base_orto(basis, rho0)
    x0 = [-scalar_prod(b, K, rho0).real for b in basis]
    res = opt.minimize(test, x0, args=(rho, basis))
    k = sum([-u*b for u,b in zip(res.x, basis)])        
    sigma = (.5*(k+k.dag())).expm()
    sigma = sigma/sigma.tr()
    return sigma
    
    
def thermal_projected_rho(rho, basis, rho0):
    # basis = base_orto(basis, rho0)
    sigma = proj_op(logM(rho), basis, rho0).expm()
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
    
  
    
def mesolve_maxent(H, rho0, tlist, c_ops=None, e_ops=None, 
                reduce_function=None,
                args=None, options=None, 
                progress_bar=None, _safe_mode=True):
    """
    Esta función evoluciona rho0 de acuerdo a la ecuación de Lindblad con 
    los parámetros indicados, pero aplica a cada tiempo en tlist 
    la función `reduce_function`. Esta funcion reduce el estado a un estado
    de la familia Max-Ent, y continúa hasta el siguiente valor de t.
    
    Por ejemplo, si se pasa como `reduce_function` la función
    ```
    def reduce_product(rho):
        res = qutip.tensor(rho.ptrace([0]),rho.ptrace([1]))
        res = res/res.tr()
        return res
    ```
    el estado evolucionará en todo tiempo como un estado producto, siguiendo una dinámica
    semejante a la de Lindblad.
    
    """
    
    # Inicialización: dejando las cosas como las quiere qutip
    if isinstance(e_ops, qutip.qobj.Qobj):
        e_ops = [e_ops]

    if isinstance(e_ops, dict):
        e_ops_dict = e_ops
        e_ops = [e for e in e_ops.values()]
    else:
        e_ops_dict = None
        
    if progress_bar is None:
        progress_bar = qutip.ui.progressbar.BaseProgressBar()
    elif progress_bar is True:
        progress_bar = qutip.ui.progressbar.TextProgressBar()
    
    
    if options is None:
        options = qutip.solver.Options()
    
    if args is None:
        args = {}

    # Start: first reduce the state and evaluate e_ops for the initial state
    progress_bar.start(len(tlist))
    
    states = [rho0] if e_ops is None else []

    if reduce_function:
        rho0 = reduce_function(rho0)

    idx_t = 0
    if e_ops:
        if type(e_ops) in (tuple, list):
            expects = [[qutip.expect(rho0, op) for op in e_ops]]
        else: # is a function?
            e_ops(tlist[0], rho0)
            expects = []
    else:
        expects = []
    
    for t1, t2 in zip(tlist[:-1],tlist[1:]):
        # Aquí resolvemos la ecuación en el intervalito
        res_inst = qutip.mesolve(H, rho0, [t1,t2], c_ops=c_ops,
                                e_ops = None, args=None, options=None,
                                progress_bar=None, _safe_mode=_safe_mode)
        rho0 = res_inst.states[-1]
        # reducimos...
        if reduce_function:
            rho0 = reduce_function(rho0)
        # y salvamos
        if e_ops:
            if type(e_ops) in (tuple, list):
                expects.append([qutip.expect(rho0, op) for op in e_ops])
            else: # is a function?
                e_ops(t2, rho0)
                
        idx_t = idx_t + 1
        progress_bar.update(idx_t)

            
    progress_bar.finished()

    if e_ops_dict:
        res_inst.expect = {e: expects[n] for n, e in enumerate(e_ops_dict.keys())}
    else:
        res_inst.expect = expects

    res_inst.states = states
    res_inst.times = tlist    
    return res_inst

# Ejemplo: Dos bosones acoplados entre sí en un sistema cerrado vs su aproximación Max Ent

#Parámetros:
omega1 = 3
omega2 = np.sqrt(48)
temp = 1
temp1 = .3
temp2 = 1.
gaussian = False
deltat = 200.

sampling = int(.25*max(1,omega1, omega2)*deltat)
ts = np.linspace(0, deltat, sampling)


distance = bures
dim1=10
dim2=10

# Base de operadores
basis_op_bos1 = [qutip.qeye(dim1), qutip.create(dim1),qutip.create(dim1).dag(),qutip.num(dim1)]
basis_op_bos2 = [qutip.qeye(dim2), qutip.create(dim2),qutip.create(dim2).dag(),qutip.num(dim2)]
prod_basis_op = prod_basis(basis_op_bos1,basis_op_bos2)
# operadores útiles
a1 = basis_op_bos1[1]
a2 = basis_op_bos2[1]
n1 = basis_op_bos1[3]
n2 = basis_op_bos2[3]
id1 = basis_op_bos1[0]
id2 = basis_op_bos2[0]

#sum_basis_op = [qutip.tensor(n1,id2), qutip.tensor(id1,n2)]
sum_basis_op = [qutip.tensor(op1,id2) for op1 in basis_op_bos1] + [qutip.tensor(id1,op2) for op2 in basis_op_bos2]

basis_op_bos_corr = prod_basis_op
"""
This basis works for U(1) invariant initial states
basis_op_bos_corr = [qutip.tensor(id1, id2),
                     qutip.tensor(id1, n2),
                     qutip.tensor(n1, id2),
                     qutip.tensor(a1, a2.dag())+qutip.tensor(a1.dag(), a2),
                     1j*(qutip.tensor(a1, a2.dag())-qutip.tensor(a1.dag(), a2)),
                    ]
"""

##########  Hamiltonianos ############33

alpha = 0.05

H_bos1 =  omega1 * qutip.tensor(n1, id2)
H_bos2 =  omega2 * qutip.tensor(id1, n2)
H_i = alpha * qutip.tensor(a1+a1.dag(), a2+a2.dag())
H = H_bos1 + H_bos2 + H_i

######   Estado inicial ################
TOp1 = (.5*(a1.dag()-a1)).expm()
TOp2 = (.5*(a2.dag()-a2)).expm()
rho0 = qutip.tensor(TOp1*(-n1/temp1).expm()*TOp1.dag(), TOp2*(-n2/temp2).expm()*TOp2.dag())
rho0 = rho0/rho0.tr()

# Base ortogonal respecto a rho0
basis_op_ortho_prod = base_orto(prod_basis_op, rho0)
basis_op_ortho_sum = base_orto(sum_basis_op, rho0)

print(datetime.datetime.now())

# Evolución exacta del estado
print("exact")
exact_states = []
qutip.mesolve(H, rho0, ts, e_ops=lambda t,rho: exact_states.append(rho))
with open("exact_states-T.pkl", "wb") as f:
    pickle.dump(exact_states, f)

print(datetime.datetime.now())


# Evolución max-ent del estado (sin correlaciones)

print("descorrelacionado")
def descorrelacionar(rho):
    sigma = qutip.tensor(rho.ptrace([0]), rho.ptrace([1]))
    sigma = sigma/sigma.tr()
    return sigma

print("tp")

def one_body_me_tp(rho):
    #prod_basis_op
    print("  step ", datetime.datetime.now())
    sigma = thermal_projected_rho(rho,basis_op_ortho_sum, rho0)
    sigma = sigma / sigma.tr()
    with open("me_1_tp_states-T.pkl", "wb") as f:
        pickle.dump(me_1_tp_states, f)
    return sigma

me_1_tp_states = []
mesolve_maxent(H, rho0, ts, e_ops=lambda t,rho: me_1_tp_states.append(rho), reduce_function=one_body_me_tp)

print("me")


def one_body_me(rho):
    #prod_basis_op
    print("  step ", datetime.datetime.now())
    sigma = maxent_rho(rho, basis_op_ortho_sum)
    # sigma = thermal_projected_rho(rho,basis_op_ortho_sum, rho0)
    sigma = sigma / sigma.tr()
    with open("me_1_full_states-T.pkl", "wb") as f:
        pickle.dump(me_1_states, f)
    return sigma

me_1_states = []
mesolve_maxent(H, rho0, ts, e_ops=lambda t,rho: me_1_states.append(rho), reduce_function=one_body_me)


# Evolución max-ent del estado (correlaciones de pares)
print("correlaciones de pares")
print("tp")
def pairwise_me_tp(rho):
    #prod_basis_op
    print("  step ", datetime.datetime.now())
    sigma = thermal_projected_rho(rho,basis_op_ortho_prod, rho0)
    sigma = sigma / sigma.tr()
    with open("me_2_tp_states-T.pkl", "wb") as f:
        pickle.dump(me_2_tp_states, f)
    return sigma

me_2_tp_states = []
mesolve_maxent(H, rho0, ts, e_ops=lambda t,rho: me_2_tp_states.append(rho), reduce_function=pairwise_me_tp)

print("me")
def pairwise_me(rho):
    #prod_basis_op
    print("  step ", datetime.datetime.now())
    sigma = maxent_rho(rho, basis_op_ortho_prod)
    #sigma = thermal_projected_rho(rho,basis_op_ortho_prod, rho0)
    sigma = sigma / sigma.tr()
    with open("me_2_full_states-T.pkl", "wb") as f:
        pickle.dump(me_2_states, f)
    return sigma

me_2_states = []
mesolve_maxent(H, rho0, ts, e_ops=lambda t,rho: me_2_states.append(rho), reduce_function=pairwise_me)


### Plots

fig = plt.figure()

ax2 = plt.subplot(2,2,2)
ax1 = plt.subplot(2,2,4, sharey=ax2)
ax3 = plt.subplot(1,2,1)

ax2.xaxis.tick_top()
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
fig.tight_layout(pad=-.5)


ax3.set_ylabel("$S_{rel,AB}$")
ax3.set_xlabel("t[s]")
ax3.plot(ts, [rel_entropy(rho, sigma)  for rho, sigma in zip(exact_states, me_1_states)], color="blue" , label="1-body MaxEnt")
ax3.plot(ts, [rel_entropy(rho, sigma)  for rho, sigma in zip(exact_states, me_2_states)], color="orange", label="2-body MaxEnt")
ax3.plot(ts, [rel_entropy(rho, sigma)  for rho, sigma in zip(exact_states, me_1_tp_states)], color="red", label="1-body TP")
ax3.plot(ts, [rel_entropy(rho, sigma)  for rho, sigma in zip(exact_states, me_2_tp_states)], color="green", label="2-body TP")
ax3.legend()


ax2.set_ylabel("$S_{rel,A}$")
ax2.plot(ts, [rel_entropy(rho.ptrace([0]), sigma.ptrace([0]))  for rho, sigma in zip(exact_states, me_1_states)], color="blue", label="me1")
ax2.plot(ts, [rel_entropy(rho.ptrace([0]), sigma.ptrace([0]))  for rho, sigma in zip(exact_states, me_2_states)], color="orange",  label="me2")
ax2.plot(ts, [rel_entropy(rho.ptrace([0]), sigma.ptrace([0]))  for rho, sigma in zip(exact_states, me_1_tp_states)], color="red", label="TP1")
ax2.plot(ts, [rel_entropy(rho.ptrace([0]), sigma.ptrace([0]))  for rho, sigma in zip(exact_states, me_2_tp_states)], color="green", label="TP2")


ax1.set_xlabel("t[s]")
ax1.set_ylabel("$S_{rel,B}$")
ax1.plot(ts, [rel_entropy(rho.ptrace([1]), sigma.ptrace([1]))  for rho, sigma in zip(exact_states, me_1_states)], color="blue", label="me1")
ax1.plot(ts, [rel_entropy(rho.ptrace([1]), sigma.ptrace([1]))  for rho, sigma in zip(exact_states, me_2_states)], color="orange", label="me2")
ax1.plot(ts, [rel_entropy(rho.ptrace([1]), sigma.ptrace([1]))  for rho, sigma in zip(exact_states, me_1_tp_states)], color="red", label="tp1")
ax1.plot(ts, [rel_entropy(rho.ptrace([1]), sigma.ptrace([1]))  for rho, sigma in zip(exact_states, me_2_tp_states)], color="green", label="tp2")

plt.savefig(f"b{dim1}xb{dim2}-RS-BXSopenNR.jpeg")


fig = plt.figure()
ax1 = plt.subplot(1,1,1)

ax1.plot(ts, [qutip.expect(rho,qutip.tensor(n1,id2))  for rho in exact_states],label="exact A")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(n1,id2))  for rho in me_1_states], color="blue", label="me 1 A")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(n1,id2))  for rho in me_2_states], color="orange",label="me 2 A")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(n1,id2))  for rho in me_1_tp_states], color="red",label="tp 1 A")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(n1,id2))  for rho in me_2_tp_states], color="green",label="tp 2 A")
ax1.legend()

plt.savefig(f"b{dim1}xb{dim2}-n1-BXSopenNR.jpeg")

fig = plt.figure()
ax1 = plt.subplot(1,1,1)

ax1.plot(ts, [qutip.expect(rho,qutip.tensor(id1,n2))  for rho in exact_states],label="exact B")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(id1,n2))  for rho in me_1_states], color="blue" , label="me 1 B")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(id1,n2))  for rho in me_2_states], color="orange", label="me 2 B")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(id1,n2))  for rho in me_1_tp_states], color="red", label="tp 1 B")
ax1.plot(ts, [qutip.expect(rho,qutip.tensor(id1,n2))  for rho in me_2_tp_states], color="green",label="tp 2 B")
ax1.legend()
plt.savefig(f"b{dim1}xb{dim2}-n2-BXSopenNR.jpeg")
