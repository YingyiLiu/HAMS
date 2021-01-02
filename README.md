<!--
  Title: HAMS
  Description: An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.
  Authors: Yingyi Liu.
  -->

# HAMS
**An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.**

[![License: Apache v2](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/body_mesh.png" width="500"></p>
  
HAMS (Hydrodynamic Analysis of Marine Structures) is an open-source software for the numerical computation of the wave effect upon marine structures. It is based on boundary integral equations in the potential flow theory for analysis of wave-structure interactions. It is currently written in FORTRAN 90. The code has been developed by the author Yingyi Liu for almost a decade. HAMS is released in the hope that it will be useful, 

HAMS is freely distributed under the Apache License, Version 2.0, http://www.apache.org/licenses/LICENSE-2.0, and may be modified and extended by researchers who intend to enhance its capabilities and port the code to other platforms. 

The success of HAMS should to a large portion be attributed to Prof. Bin Teng (Dalian University of Technology), who has tutored me the theory of potential flow in marine hydrodynamics and the programming skills using the [Boundary Element Method](https://en.wikipedia.org/wiki/Boundary_element_method). The code structure and the coding style of HAMS are exactly two of the examples that I have learned and inherited from Prof. Bin Teng.

## Theoretical Basis

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/global_coordinate_system.png" width="500"></p>
  
### - Please refer to the following papers for the theory:

* Yingyi Liu (2019). HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application. _Journal of Marine Science and Engineering_ 7: 81.

* Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He (2018). A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions. _Energy Conversion and Management_, 174: 516-536.

* Yingyi Liu, Changhong Hu, Makoto Sueyoshi, Hidetsugu Iwashita, Masashi Kashiwagi (2016). Motion response prediction by hybrid panel-stick models for a semi-submersible with bracings. _Journal of Marine Science and Technology_, 21:742–757.

Please cite the above papers in your relevant publications, reports, etc., if the HAMS code or its executable program has contributed to your work.

## Generated numerical results

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

* HAMS can import meshes containing triangular panel type, quadrilateral panel type, or both.
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/panel_local_coordinates.png" width="650"></p>

### - OpenMP parallel processing

* HAMS can be run in the parallel mode using OpenMP techniques on PC's with multiple processors (CPU's).
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/OpenMP_parallel_process.png" width="700"></p>
  
### - Computational efficiency

* HAMS can be run in the parallel mode using OpenMP techniques on PC's with multiple processors (CPU's).
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/computational_efficiency.png" width="500"></p>

## Useful Links

The following open-source software can be used to view the HAMS results: <br/>
[1]. [BEMRosetta](https://github.com/izabala123/BEMRosetta).<br/> Developed by Iñaki Zabala, Markel Peñalba, Yerai Peña-Sanchez.
[2]. [BEMIO](https://wec-sim.github.io/bemio/).<br/> Developed by National Renewable Energy Laboratory and Sandia 
Sandia National Laboratories. 


## License

HAMS is free software: you can redistribute it and/or modify it under the terms of the Apache License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.\

HAMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache License for details. You should have received a copy of the Apache License along with HAMS. If not, see http://www.apache.org/licenses/LICENSE-2.0 <br/>


