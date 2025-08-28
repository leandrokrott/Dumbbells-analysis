import numpy as np
import os
import re

# >>> DEFINIÇÕES MANUAIS DE PARÂMETROS <<<
k_fixed = "100.0"   # valor de k a ser fixado

# Função para corrigir deslocamentos considerando as condições periódicas de contorno
def unwrap_positions(positions, box_size):
    unwrapped = np.copy(positions)
    for i in range(1, len(positions)):
        delta = positions[i] - positions[i - 1]
        delta -= np.round(delta / box_size) * box_size
        unwrapped[i] = unwrapped[i - 1] + delta
    return unwrapped

# Função para calcular o MSD
def compute_msd(trajectories, timesteps):
    num_particles, num_steps, _ = trajectories.shape
    msd = np.zeros(num_steps)
    for t in range(num_steps):
        displacements = trajectories[:, t, :] - trajectories[:, 0, :]
        squared_displacements = np.sum(displacements ** 2, axis=1)
        msd[t] = np.mean(squared_displacements)
    return msd

# Lista de arquivos no diretório atual que correspondem aos parâmetros fixos
input_files = [
    f for f in os.listdir()
    if re.match(rf"dump_cm_k_{k_fixed}_r0_([\d\.]+)_L_([\d\.]+)_T_([\d\.]+)\.lammpstrj", f)
]


for input_file_name in input_files:
    print(f"Processando {input_file_name}...")

    # Extraindo parâmetros do nome do arquivo
    match = re.search(r"dump_cm_k_([\d\.]+)_r0_([\d\.]+)_L_([\d\.]+)_T_([\d\.]+)\.lammpstrj", input_file_name)
    if not match:
        print(f"Nome de arquivo {input_file_name} não está no formato esperado. Pulando.")
        continue
    k, r0, L, T = match.groups()
    box_bounds = [float(L) * 2] * 3  # como o intervalo vai de -L a L
    box_size = np.array(box_bounds)

    with open(input_file_name, "r") as file:
        lines = file.readlines()

    timesteps = []
    positions = {}
    current_timestep = None
    num_atoms = 0

    for i, line in enumerate(lines):
        if "ITEM: TIMESTEP" in line:
            current_timestep = int(lines[i + 1].strip())
            timesteps.append(current_timestep)
        elif "ITEM: NUMBER OF ATOMS" in line:
            num_atoms = int(lines[i + 1].strip())
        elif "ITEM: ATOMS" in line:
            start_idx = i + 1
            timestep_positions = np.zeros((num_atoms, 3))
            for j in range(num_atoms):
                data = lines[start_idx + j].split()
                atom_id = int(data[0]) - 1  # ajusta índice
                timestep_positions[atom_id] = np.array(list(map(float, data[2:5])))  # x, y, z
            positions[current_timestep] = timestep_positions

    timesteps.sort()
    num_steps = len(timesteps)

    trajectories = np.zeros((num_atoms, num_steps, 3))
    for t_idx, t in enumerate(timesteps):
        trajectories[:, t_idx, :] = positions[t]

    for i in range(num_atoms):
        trajectories[i] = unwrap_positions(trajectories[i], box_size)

    msd = compute_msd(trajectories, timesteps)

    output_file_msd = f"msd_cm_k_{k}_r0_{r0}_L_{L}_T_{T}.dat"
    np.savetxt(output_file_msd, np.column_stack((timesteps, msd)), comments="")
    print(f"MSD salvo em {output_file_msd}")

print("Processamento completo.")
