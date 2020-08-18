"""
In questo file ci si propone di verificare i conti sul controesempio:
conditional non-signalling =/= oCP.

I conti dell'articolo tornano a meno della formula B2.
A me torna che il vettore dello stato è: [(a-2c)/3, b/3, (2a-c)/3]

Convenzione: S x E_s x E_r
"""

tollerated_error = 1.e-12
number_of_test = 1000

from qutip import *
from qutip.qip.operations import swap
import numpy as np
from math import sqrt
import random

def vector2qubit(v):
    return 1./2.*(identity(2) + v[0]*sigmax() + v[1]*sigmay() + v[2]*sigmaz())

def qubit2vector(q):
    return np.array([(q*sigmax()).tr(), (q*sigmay()).tr(), (q*sigmaz()).tr()])

# Evoluzioni
U_sr = swap(3,[0,2])
U_ts = 1./sqrt(3.) * tensor(tensor(identity(2),identity(2)) + 1j * tensor(sigmay(),identity(2)) + 1j * tensor(sigmax(), sigmaz()), identity(2))
U_tr = U_ts * U_sr

def verify_simulation(v_r, v_s):
    # Stato iniziale del sistema
    rho_r = vector2qubit(v_r)
    eta_r = (tensor(sigmaz(),sigmax())+tensor(qeye(2),qeye(2)))/4
    system_r = tensor(rho_r, eta_r)

    # Evoluzione interrotta
    system_s_minus = U_sr * system_r * U_sr.dag()
    rho_s_minus = system_s_minus.ptrace(0)
    rho_s_plus = vector2qubit(v_s)
    eta_s = system_s_minus.ptrace([1,2])
    system_s_plus = tensor(rho_s_plus, eta_s)
    system_t_int = U_ts * system_s_plus * U_ts.dag()
    rho_t_int = system_t_int.ptrace(0)


    # Evoluzione completa
    system_t = U_tr * system_r * U_tr.dag()
    rho_t = system_t.ptrace(0)

    # Verifiche
    assert(rho_s_minus == identity(2)/2.)
    assert(rho_t == identity(2)/2. + sigmay()/3.)
    assert(eta_s == tensor(identity(2)/2., rho_r))
    v_t_theo = np.array([(v_s[0]-2*v_s[2])/3., v_s[1]/3., (2*v_s[0]-v_s[2])/3.])
    v_t = qubit2vector(rho_t_int)
    eps = np.abs(v_t_theo - v_t)
    assert(np.all(eps <= tollerated_error))

def random_point_in_sphere():
    while True:
        x = random.uniform(-1., 1.)
        y = random.uniform(-1., 1.)
        z = random.uniform(-1., 1.)
        if x**2 + y**2 + z**2 <=1.:
            return np.array([x,y,z])


for i in range(1,1000):
    print(i)
    verify_simulation(random_point_in_sphere(),random_point_in_sphere())
