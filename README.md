# Dumbbells Simulation in NVT Ensemble

This repository contains scripts and workflows for simulating dumbbells in the NVT ensemble, computing orientational and translational order parameters, excess entropy, self-diffusion, and average of monomer-monomer distance.

---

## 1. Generate Initial Configuration

Use the script `in.dumbbells_config` with the input files `alan.table` and `mymol.txt` for LAMMPS. Specify:

- Spring constant: `k_el`
- Half-box size: `L`
- Desired temperature (start with `T = 1.00`)

The initial configuration file will be generated as:

config_k_${k_el}_L${L}.dat


## 2. Run the Main Simulation

Run the main program using LAMMPS: 

in.dumbbells_loop_NVT

For each specific value of `k_el` and `L`, loop over temperatures from `T = 0.900` down to `T = 0.050`.  

Outputs:

- **Log file**: Contains simulation details
- **Dump file**: Contains trajectories and velocities of all particles

---

## 3. Compute Orientational Order Parameters

Use the script:

qs_parameters.py


- Calculates `q_l` for all `l`, focusing on `q_4` and `q_6`.

---

## 4. Radial Distribution Function (RDF) and Center of Mass Analysis

- Run `gr_CM.f90` to compute the RDF of the center of mass (and optionally all particles and/or even (or odd) particles).  
- Generate `dump_cm` files for the center-of-mass trajectories.

gfortran gr_CM.f90 -o call (to generate the executable call)

./script_gr_dumpcm (it will use the executable call and run for all variables)

---

## 5. Translational Order Parameter and Excess Entropy

- Using the RDF from `gr_CM.f90`, calculate:
  - Translational order parameter `τ`
  - Two-body excess entropy `s_2`  

Using the program:

gfortran tau_s2.f90 -o tr (to generate the executable tr)

./script_tau_s2 (it will use the executable tr and run for all variables)


---

## 6. Mean Square Displacement and Self-Diffusion

- Compute mean square displacement (MSD) of the center of mass:

MSD.py


- Compute self-diffusion coefficient:

self_diffusion.py


---

## 7. Average of monomer-monomer distance

- Compute the average of monomer-monomer distance:

lambda_average.py

---

## References

- Scripts and inputs: `alan.table`, `mymol.txt`
- Main programs: `in.dumbbells_config`, `in.dumbbells_loop_NVT`
- Analysis: `qs_parameters.py`, `gr_CM.f90`, `tau_s2.f90`, `MSD.py`, `self_diffusion.py`, `lambda_average.py`




