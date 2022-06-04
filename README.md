# GSD_BAOA_BAOAB

this repository contains supporting information to our publication:<br>
S.Kieninger, B.G. Keller, **"GROMACS Stochastic Dynamics and BAOAB are equivalent configurational sampling algorithms"**, arXiv preprint:  	arXiv:2204.02105 (2022) <br>
<br>
In our publication, we studied a one-dimensional model system, an ideal gas and TIP3P bulk water with OpenMM<sup>[[1]](#1)[[2]](#2)</sup> (version 7.4.2) in combination with OpenMMTools<sup>[[3]](#3)</sup>, GROMACS<sup>[[4]](#4)</sup> (version 2019.4) and Python3. 
You can find the corresponding input files, programme codes and simulation protocols of the these systems in this repository. 
Additionally, we provide all results presented in the figures as files in npy-format (```numpy.lib.format```).


<br>
<br>
<h2>Overview </h2>

<h4>1D model system </h4>

All MD simulations and the post analysis were performed with Python3. 
We included Python scripts on how to load the provided results.
The directory tree is as follows:

```
1D_model_system
  |
  ├── distributions 
  |       |   compute_distributions.py
  |       |   compute_relative_error_temperatures.py
  |       |   compute_temperatures.py
  |       |
  |       └── results
  |              how_to_load_distributions.py
  |              how_to_load_relative_error_temperatures.py
  |              results as .npy
  |              logs as .txt
  |
  └── example_paths
          |   generate_example_paths.py
          |   random_numbers_example_paths.npy
          |
          └── results
                 how_to_load_results.py
                 results as .npy
```

<h4>Ideal gas </h4>

All MD simulations were performed with OpenMM and the post analysis was carried out with Python3.
We included Python scripts on how to load the provided results.
The directory tree is as follows:

```
ideal_gas
  |
  └── thermal_rates
          |
          ├── OpenMM
          |       mdrun_ABOBA.py
          |       mdrun_GSD.py
          |
          └── results
                  how_to_load_thermal_rates.py
                  results as .npy
```

<h4>TIP3P bulk water </h4>

We performed MD Simulations with OpenMM as well as GROMACS. The post analysis was carried out with Python3.
We included Python scripts on how to load the provided results.
The directory tree is as follows:

```
water_box
  |
  └── energy_distributions
          |
          ├── structure_and_toplogy
          |       water_box.gro
          |       water_box.top
          |
          ├── Gomacs
          |       simulation.mdp
          |
          ├── OpenMM
          |       mdrun_ABOBA.py
          |       mdrun_GSD.py
          |
          └── results
                  how_to_load_distributions.py
                  results as .npy
```
<br>
<h2>References </h2>

<a id="1">[1]</a>
P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp, L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora, B. R. Brooks, V. S. Pande,
*PLoS Comput. Biol.*, **13**, (2017), 1 

<a id="2">[2]</a>
<https://github.com/openmm/openmm>

<a id="3">[3]</a>
<https://github.com/choderalab/openmmtools>


<a id="4">[4]</a>
H. Berendsen, D. van der Spoel, R. van Drunen,
*Comput. Phys. Commun.*, **91**, (1995), 43-65 
