<!--
  Title: HAMS
  Description: An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.
  Authors: Yingyi Liu.
  -->

# HAMS
**An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.**

[![License: Apache v2](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

HAMS (Hydrodynamic Analysis of Marine Structures) is an open-source software for the numerical computation of the wave effect upon marine structures. It is based on boundary integral equations in the potential flow theory for analysis of wave-structure interactions. It is currently written in FORTRAN 90. The code has been developed by the author Yingyi Liu for around a decade. 

HAMS is freely distributed under the Apache License, Version 2.0, http://www.apache.org/licenses/LICENSE-2.0, and may be modified and extended by researchers who intend to enhance its capabilities and port the code to other platforms.

It should be noted that, any modified version should be licensed under the LGPL License and be released open-publicly as well. The contributors can add their names in the "contributors list" ahead of the modified subroutine(s).

[Boundary Element Methods](https://en.wikipedia.org/wiki/Boundary_element_method) are extensively used to model hydrodynamic forces in offshore devices like ships, offshore wind platforms and wave energy converters. These solvers use device geometry mesh to get some hydrodynamics coefficients as radiation damping, added mass, wave diffraction force, and wave excitation force. All these data is saved in file formats incompatible between them. These may avoid to use the coefficients between programs. 

BEMRosetta allows to load the hydrodynamic coefficients from a format saving it in another. In addition it allows to compare the results obtained between programs, the results between similar geometries and the same geometry with different discretization levels.

Moreover, BEMRosetta allows to view and visually compare the meshes from different programs.

BEMRosetta runs on Windows and Linux, **no install is necessary in Windows** [(see Install)](https://github.com/izabala123/BEMRosetta/tree/master/install), and it includes a GUI and [a command line version](https://github.com/izabala123/BEMRosetta/blob/master/other/test). 

## Theoretical Basis

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/global_coordinate_system.png" width="500" title="Free-surface elevation"></p>
  
### - Please refer to the following papers:

Yingyi Liu (2019). HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application. Journal of Marine Science and Engineering 7: 81.

Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He (2018). A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions. Energy Conversion and Management, 174: 516-536.

Please cite the above papers in your relevant publications if the HAMS code or its executable program has contributed to your work.

## Generated results

### - Hydrodynamic coefficients

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/hydrodynamic_coefficients.png" width="900"></p>

### - Wave excitation force

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/wave_excitation_force_plot.png" width="700"></p>

### - Motion RAOs

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/RAO_of_motion.png" width="600"></p>

### - Free-surface elevation

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/free_surface_elevation.png" width="500"></p>

## Features

### - Mesh element type

  Can be triangular panel, quadrilateral panel, or both.
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/panel_local_coordinates.png" width="600"></p>

### - OpenMP parallel processing

  HAMS can be run in parallel mode on PC's with multiple processors (CPU's).
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/OpenMP_parallel_process.png" width="900"></p>
### - Computational efficiency

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/computational_efficiency.png" width="600"></p>

## Useful Links

The following open-source software can be used to view the HAMS results: <br/>
[1].[BEMRosetta](https://github.com/izabala123/BEMRosetta).<br/>
[2].[BEMIO](https://wec-sim.github.io/bemio/).<br/>


## License

HAMS is free software: you can redistribute it and/or modify it under the terms of the Apache License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.\

HAMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>


