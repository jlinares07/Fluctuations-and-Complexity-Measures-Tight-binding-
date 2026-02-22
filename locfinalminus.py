import qutip as qt
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from tqdm import tqdm

# Parámetros
N = 50
a = 1.0
hbar = 1.0
lambda_ = 0.01
U = -100.0 * lambda_
F = 4 * lambda_
V0 = 0.0 
E0 = -0.01
delta_t = 0.02
total_time = 1000.0
M=3

# Dimensión de la base
dim = (2 * N + 1) ** 2 
states = [(l, m) for l in range(-N, N + 1) for m in range(-N, N + 1)]
state_index = {state: i for i, state in enumerate(states)}

# Matriz Sparse que facilita el trabajo con N grande
H_data = sp.lil_matrix((dim, dim), dtype=complex)

for l, m in states:
    i = state_index[(l, m)]
    if l == m:
        coulomb_interac = U
    else:
        coulomb_interac = V0 / abs(l - m) if abs(l-m) != 0 else 0
    
    E_lm = E0 + coulomb_interac + F * (l + m) / hbar
    H_data[i, i] = E_lm

for l, m in states:
    i = state_index[(l, m)]
    if l < N:
        lp = l + 1
        j = state_index.get((lp, m))
        if j is not None:
            H_data[i, j] = -lambda_
            H_data[j, i] = -lambda_
    if m < N:
        mp = m + 1
        j = state_index.get((l, mp))
        if j is not None:
            H_data[i, j] = -lambda_
            H_data[j, i] = -lambda_

H_data = H_data.tocsr() #Conversión a CSR para eficiencia

##### Funciones
def calculate_system_entropy_desequilibrium_complexity(probabilities):
    prob_nonzero = probabilities[probabilities > 0]
    entropy = -np.sum(prob_nonzero * np.log(prob_nonzero)) if prob_nonzero.size > 0 else 0
    lin_entropy = 1 - np.sum(prob_nonzero**2)
    equiprobability = 1 / dim
    desequilibrium = np.sum((probabilities - equiprobability) ** 2)
    complexity = entropy * desequilibrium
    return entropy, desequilibrium, complexity, lin_entropy

def probability_current(psi_nm, states, state_index, lambda_, a, hbar):
    current = 0.0
    for l, m in states:
        i = state_index[(l, m)]
        psi_lm = psi_nm[i]
        psi_lm_minus_l = 0.0 if (l - 1, m) not in state_index else psi_nm[state_index[(l - 1, m)]].conj()
        psi_lm_minus_m = 0.0 if (l, m - 1) not in state_index else psi_nm[state_index[(l, m - 1)]].conj()
        product = psi_lm * (psi_lm_minus_l + psi_lm_minus_m)
        current += product.imag * lambda_ * 2 * (a / hbar)
    return -current

def centro_masa(probabilities, states):
    centro = 0.0
    for state_index, (l,m) in enumerate(states):
        centro += probabilities[state_index]*(l+m) / 2
    return centro
#####

# Convertir a Qobj para trabajar con sesolve
H = qt.Qobj(H_data, dims=[[dim], [dim]])
#print(H.isherm, "fase inicial")
t_list = np.linspace(0, total_time, int(total_time/delta_t))
options = {'store_states': True, 'progress_bar': 'enhanced'}
#psi0_l = qt.basis(dim, 1)
#psi0_m = qt.basis(dim, 2) #Posibilidad de trabajar con la funcion tensor, sin embargo aparentemente dificulta el trabajo computacional
#psi0 = qt.basis(dim, state_index[-3,3]) #Esto funciona cuando se quiere un estado completamente localizado
psi0 = qt.Qobj(np.zeros((dim, 1)), dims=[[dim], [1]])
psi0_loc = (qt.basis(dim, state_index[1,-1]) + qt.basis(dim, state_index[-1,1])) / np.sqrt(2)
num_states = 0
for l in range(-M, M + 1):
    for m in range(-M, M + 1):
        if (l, m) in state_index:
            idx = state_index[(l, m)]
            psi0 += qt.basis(dim, idx)
            num_states += 1
psi0 = psi0 / np.sqrt(num_states)
# Tanto psi0 como H disponen de la misma dimensionalidad, si hace falta corroborar se puede revisar rápidamente
#print(H, psi0) para revisar dimensiones
evolution = qt.sesolve(H, psi0_loc, t_list, options = options)
#H_dense = H.full()  # Usar solo para N pequeño o análisis específico de dimensiones
steps = int(total_time/delta_t)
entropies, desequilibria, complexities, currents, cm, linear_entropy, complex_sdl = [], [], [], [], [], [], []
#print(H.isherm, "fase 1")

for s in tqdm(range(steps), desc="Calculating quantities"):
    probabilities = np.array([abs(evolution.states[s].overlap(qt.basis(dim,i)))**2 for i in range(dim)]) #Agiliza la evolución de la función de onda
    psi_nm = evolution.states[s].full().flatten()
    entropy, desequilibrium, complexity, lin_entropy = calculate_system_entropy_desequilibrium_complexity(probabilities)
    current = probability_current(psi_nm, states, state_index, lambda_, a, hbar)
    centro = centro_masa(probabilities, states)
    entropies.append(entropy)
    linear_entropy.append(lin_entropy)
    desequilibria.append(desequilibrium)
    complexities.append(complexity)
    cm.append(centro)
    currents.append(current)

smax = np.max(linear_entropy)
for s in range(steps):
    delocalization = linear_entropy[s] / smax
    complexity_sdl = delocalization * (1 - delocalization)
    complex_sdl.append(complexity_sdl)


np.savez(f"Localized,U={U},V0={V0}.npz", t_list=t_list, entropies=entropies, desequilibria=desequilibria, complexities=complexities, currents = currents, cm=cm, linear_entropy=linear_entropy, complex_sdl=complex_sdl)
