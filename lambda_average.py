import numpy as np
import re
import glob
import os

def extrair_parametros(nome_arquivo):
    padrao = r"dump_k_([0-9eE.+-]+)_r0_([0-9eE.+-]+)_L_([0-9eE.+-]+)_T_([0-9eE.+-]+)\.lammpstrj"
    m = re.search(padrao, nome_arquivo)
    if m:
        return m.groups()
    else:
        raise ValueError("Nome do arquivo não segue o padrão esperado.")

def calcular_distancias_por_timestep(nome_arquivo):
    distancias_por_dimero = []
    n_dimeros = None

    with open(nome_arquivo, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break

            if line.startswith("ITEM: TIMESTEP"):
                f.readline()  # pula o timestep

                # Número de átomos
                line = f.readline()
                if not line.startswith("ITEM: NUMBER OF ATOMS"):
                    raise ValueError("Formato inesperado: esperado 'ITEM: NUMBER OF ATOMS'")
                n_atoms = int(f.readline().strip())

                if n_atoms % 2 != 0:
                    raise ValueError("Número de átomos não é par (não forma dímeros).")
                n_dimeros = n_atoms // 2

                # Caixa
                line = f.readline()
                if not line.startswith("ITEM: BOX BOUNDS"):
                    raise ValueError("Formato inesperado: esperado 'ITEM: BOX BOUNDS'")
                bounds_x = f.readline().split()
                bounds_y = f.readline().split()
                bounds_z = f.readline().split()
                Lx = float(bounds_x[1]) - float(bounds_x[0])
                Ly = float(bounds_y[1]) - float(bounds_y[0])
                Lz = float(bounds_z[1]) - float(bounds_z[0])

                # Coords
                line = f.readline()
                if not line.startswith("ITEM: ATOMS"):
                    raise ValueError("Formato inesperado: esperado 'ITEM: ATOMS...'")

                coords = []
                for _ in range(n_atoms):
                    parts = f.readline().split()
                    x, y, z = map(float, parts[2:5])
                    coords.append([x, y, z])
                coords = np.array(coords)

                delta = coords[0::2] - coords[1::2]
                delta[:, 0] -= Lx * np.round(delta[:, 0] / Lx)
                delta[:, 1] -= Ly * np.round(delta[:, 1] / Ly)
                delta[:, 2] -= Lz * np.round(delta[:, 2] / Lz)

                distancias = np.linalg.norm(delta, axis=1)

                if not distancias_por_dimero:
                    distancias_por_dimero = [[] for _ in range(n_dimeros)]

                for i in range(n_dimeros):
                    distancias_por_dimero[i].append(distancias[i])

    medias = [np.mean(d) for d in distancias_por_dimero]
    erros = [np.std(d, ddof=1) / np.sqrt(len(d)) for d in distancias_por_dimero]

    return medias, erros, n_dimeros

def processar_arquivo(nome_arquivo):
    try:
        k, r0, L, T = extrair_parametros(nome_arquivo)
        nome_saida = f"separacao_media_check_k_{k}_r0_{r0}_L_{L}_T_{T}.dat"
        checkpoint_file = f"{nome_saida}.done"

        # Se o arquivo .done existir, pular
        if os.path.exists(checkpoint_file):
            print(f"[→] Ignorado (checkpoint existe): {nome_arquivo}")
            return

        medias, erros, n_dimeros = calcular_distancias_por_timestep(nome_arquivo)

        with open(nome_saida, "w") as f:
            for i, (m, e) in enumerate(zip(medias, erros)):
                f.write(f"{i} {m:.6f} {e:.6f}\n")

        # Cria o checkpoint
        with open(checkpoint_file, "w") as f:
            f.write("OK\n")

        print(f"[✔] Arquivo processado: {nome_saida}")
    except Exception as e:
        print(f"[✘] Erro ao processar {nome_arquivo}: {str(e)}")

# Processa todos os arquivos que seguem o padrão
for arquivo in glob.glob("dump_k_*.lammpstrj"):
    processar_arquivo(arquivo)
