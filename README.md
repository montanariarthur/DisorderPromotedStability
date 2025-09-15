# Disorder-promoted stability in complex systems
Code for the analysis of disorder in several second-order dynamical systems, including Kuramoto oscillators, metamaterials, phase-amplitude oscillators, and multi-agent systems.

See the references below for more details. If you use this code, please cite:
  - AN Montanari, P Zanin, AE Motter. Disorder-promoted stability. (2025)

# Usage

- `main_dps_dynsystems.m` : Calculates the stability impromevement promoted by an *optimal* parameter perturbation for a wide variety of dynamical systems (see file for the list of systems). The optimal parameter perturbation is calculated via numerical optimization, assuming both constrained (homogeneous) and unconstrained (heterogeneous) formulations.

- `main_dps_networks.py` : Calculates the stability improvement promoted by *random* parameter perturbations over an ensemble of random network models (small-world, Erdos-Renyi, scale-free) as well as real-world network datasets. The adjacency matrix can be readily specified by the user in the code.

- `main_dps_ecology.m` : Calculates the stability improvement by *random* perturbations to the network structure on explicit (Lotka-Volterra) and implicit (linear) models.

- `main_arnoldtongues.m` : Calculates the Arnold tongues (i.e., stability regions) for pairs of coupled oscillators, considering the following models: first-order Kuramoto, first-order leaky integrator, second-order Kuramoto, and phase-amplitude oscillator.

# Dependences

- `Model` : This folder contains codes to generate the network structures, calculate the Lyapunov exponents, and numerically find the equilibrium points of dynamical systems.

All codes were tested and run in MATLAB 2023a. To run the codes, download all files in this repository to a folder of your choice and run one of the `main` scripts of your choice. All codes generate/include the required data to run the simulations and optimization; simulations can take a few minutes on a standard laptop.

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".
