import os
import numpy as np
import freud as fd

# Pasta atual
output_folder = 'output_media'
os.makedirs(output_folder, exist_ok=True)

checkpoint_file = os.path.join(output_folder, 'checkpoint.txt')

def get_dimer_com(positions, box_size):
    N = len(positions)
    if N % 2 != 0:
        raise ValueError("Número ímpar de partículas")
    coms = []
    Lx, Ly, Lz = box_size
    box = fd.box.Box(Lx, Ly, Lz)
    for i in range(0, N, 2):
        pos1 = positions[i]
        pos2 = positions[i+1]
        delta = box.wrap(pos2 - pos1)
        com = box.wrap(pos1 + delta/2)
        coms.append(com)
    return np.array(coms)

def read_timesteps(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    timesteps = []
    box_sizes = []
    i = 0
    while i < len(lines):
        if "ITEM: TIMESTEP" in lines[i]:
            i += 1
            timestep = int(lines[i].strip())
            i += 1
            if "ITEM: NUMBER OF ATOMS" not in lines[i]:
                raise ValueError("Formato inesperado")
            i += 1
            num_atoms = int(lines[i].strip())
            i += 1
            if "ITEM: BOX BOUNDS" not in lines[i]:
                raise ValueError("Formato de caixa inesperado")
            i += 1
            xlo, xhi = map(float, lines[i].strip().split())
            i += 1
            ylo, yhi = map(float, lines[i].strip().split())
            i += 1
            zlo, zhi = map(float, lines[i].strip().split())
            i += 1
            box_size = np.array([xhi - xlo, yhi - ylo, zhi - zlo])
            box_sizes.append(box_size)
            if "ITEM: ATOMS" not in lines[i]:
                raise ValueError("Formato inesperado")
            i += 1
            data = []
            for _ in range(num_atoms):
                parts = lines[i].strip().split()
                data.append([float(p) for p in parts])
                i += 1
            timesteps.append(np.array(data))
        else:
            i += 1
    return timesteps, box_sizes

def process_single_timestep(positions, box_size):
    Lx, Ly, Lz = box_size
    box = fd.box.Box(Lx, Ly, Lz)
    vor = fd.locality.Voronoi()
    nlist = vor.compute((box, positions)).nlist
    q_particles = np.array([fd.order.Steinhardt(l=l, average=True).compute(
        (box, positions), nlist).particle_order for l in range(3, 13)])
    even_pos = positions[::2]
    nlist_even = vor.compute((box, even_pos)).nlist
    q_even = np.array([fd.order.Steinhardt(l=l, average=True).compute(
        (box, even_pos), nlist_even).particle_order for l in range(3, 13)])
    com_pos = get_dimer_com(positions, box_size)
    nlist_com = vor.compute((box, com_pos)).nlist
    q_com = np.array([fd.order.Steinhardt(l=l, average=True).compute(
        (box, com_pos), nlist_com).particle_order for l in range(3, 13)])
    return {'particles': q_particles.T, 'even': q_even.T, 'com': q_com.T}

def process_file(file_path):
    timestep_data, box_sizes = read_timesteps(file_path)
    filename = os.path.splitext(os.path.basename(file_path))[0]
    all_particles, all_even, all_com = [], [], []
    for data, box_size in zip(timestep_data, box_sizes):
        positions = data[:, 2:5]
        results = process_single_timestep(positions, box_size)
        all_particles.append(results['particles'])
        all_even.append(results['even'])
        all_com.append(results['com'])
    avg_particles = np.mean(np.array(all_particles), axis=0)
    avg_even = np.mean(np.array(all_even), axis=0)
    avg_com = np.mean(np.array(all_com), axis=0)
    np.savetxt(os.path.join(output_folder, f'q-particles-M-{filename}.dat'),
               avg_particles, fmt='%14.12f')
    np.savetxt(os.path.join(output_folder, f'q-even-M-{filename}.dat'),
               avg_even, fmt='%14.12f')
    np.savetxt(os.path.join(output_folder, f'q-com-M-{filename}.dat'),
               avg_com, fmt='%14.12f')

# Checkpoint
if os.path.exists(checkpoint_file):
    with open(checkpoint_file, 'r') as f:
        processed_files = set(line.strip() for line in f if line.strip())
else:
    processed_files = set()

print("Processando arquivos na pasta atual...")
for file in sorted(os.listdir('.')):
    if file.startswith('dump_k_') and '_r0_1.0_' in file and file.endswith('.lammpstrj'):
        if file in processed_files:
            print(f"  Pulando {file} (já processado)")
            continue
        print(f"\nArquivo: {file}")
        try:
            process_file(file)
            print("  Processado com sucesso")
            with open(checkpoint_file, 'a') as f:
                f.write(f"{file}\n")
        except Exception as e:
            print(f"  ERRO: {e}")
            with open(os.path.join(output_folder, 'error_log.txt'), 'a') as log:
                log.write(f"\nErro em {file}:\n{e}\n")

print("\nProcessamento concluído!")
