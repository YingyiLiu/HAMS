<!--
  Title: HAMS
  Description: An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.
  Authors: Yingyi Liu.
  -->

# HAMS
**An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.**

[![License: Apache v2](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/body_mesh.png" width="500"></p>
  
HAMS (Hydrodynamic Analysis of Marine Structures) is an open-source software for the numerical computation of the wave effect upon marine structures. It is based on boundary integral equations in the potential flow theory for analysis of wave-structure interactions. It is currently written in FORTRAN 90. The code has been developed by the author Yingyi Liu for around a decade. 

HAMS is freely distributed under the Apache License, Version 2.0, http://www.apache.org/licenses/LICENSE-2.0, and may be modified and extended by researchers who intend to enhance its capabilities and port the code to other platforms.

It should be noted that, any modified version should be licensed under the LGPL License and be released open-publicly as well. The contributors can add their names in the "contributors list" ahead of the modified subroutine(s).

## Theoretical Basis

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/global_coordinate_system.png" width="500"></p>
  
### - Please refer to the following papers for the theory:

* Yingyi Liu (2019). HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application. Journal of Marine Science and Engineering 7: 81.

* Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He (2018). A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions. Energy Conversion and Management, 174: 516-536.

Please cite the above papers in your relevant publications if the HAMS code or its executable program has contributed to your work.

## Generated results

### - Hydrodynamic coefficients

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/hydrodynamic_coefficients.png" width="800"></p>

### - Wave excitation force

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/wave_excitation_force_plot.png" width="550"></p>

### - Motion RAOs

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/RAO_of_motion.png" width="500"></p>

### - Free-surface elevation

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/free_surface_elevation.png" width="450"></p>

## Features

### - Mesh element type

* Can be triangular panel, quadrilateral panel, or both.
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/panel_local_coordinates.png" width="600"></p>

### - OpenMP parallel processing

* HAMS can be run in parallel mode on PC's with multiple processors (CPU's).
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/OpenMP_parallel_process.png" width="700"></p>
  
### - Computational efficiency

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/computational_efficiency.png" width="500"></p>

## Useful Links

The following open-source software can be used to view the HAMS results: <br/>
[1]. [BEMRosetta](https://github.com/izabala123/BEMRosetta).<br/>
[2]. [BEMIO](https://wec-sim.github.io/bemio/).<br/>


## License

HAMS is free software: you can redistribute it and/or modify it under the terms of the Apache License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.\

HAMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>


