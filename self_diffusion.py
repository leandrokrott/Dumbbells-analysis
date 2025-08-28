import numpy as np
import glob
import re
import os
from collections import defaultdict

def calcular_difusao(time, r2):
    coeffs = np.polyfit(time, r2, 1)  # regressão linear
    dr2_dt = coeffs[0]
    return dr2_dt / 6

def extrair_parametros(nome_arquivo):
    padrao = r'msd_cm_k_([\d.]+)_r0_([\d.]+)_L_([\d.]+)_T_([\d.]+)\.dat'
    match = re.match(padrao, nome_arquivo)
    if match:
        return tuple(map(float, match.groups()))
    return None

# Dicionários para T vs D e para ρ vs D
dados_diff_k_rho = defaultdict(list)  # chave: (k, rho) -> lista (T, D)
dados_diff_k_T = defaultdict(list)    # chave: (k, T) -> lista (rho, D)

for arquivo in glob.glob("msd_cm_k_*.dat"):
    parametros = extrair_parametros(os.path.basename(arquivo))
    if not parametros:
        continue

    k, r0, L, T = parametros
    rho = 2000.0 / (2 * L) ** 3

    dados = np.loadtxt(arquivo)
    tempo, r2 = dados[:, 0], dados[:, 1]
    D = calcular_difusao(tempo, r2)

    dados_diff_k_rho[(k, rho)].append((T, D))
    dados_diff_k_T[(k, T)].append((rho, D))

# Salvar T vs D (já existente)
for (k, rho), lista in dados_diff_k_rho.items():
    lista.sort(key=lambda x: x[0])  # ordenar por T
    with open(f"T_vs_D_k_{k}_rho_{rho:.4f}.dat", 'w') as f:
        for T, D in lista:
            f.write(f"{T:.5f} {D:.8f}\n")

# Salvar ρ vs D para cada T
for (k, T), lista in dados_diff_k_T.items():
    lista.sort(key=lambda x: x[0])  # ordenar por rho
    with open(f"RHO_vs_D_k_{k}_T_{T:.3f}.dat", 'w') as f:
        for rho, D in lista:
            f.write(f"{rho:.5f} {D:.8f}\n")
