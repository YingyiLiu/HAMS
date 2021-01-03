<!--
  Title: HAMS
  Description: An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.
  Authors: Yingyi Liu.
  -->

# HAMS
**An open-source computer program for the analysis of wave diffraction and radiation of three-dimensional floating or submerged structures.**

[![License: Apache v2](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/body_mesh.png" width="70%" /></p>
  
HAMS (Hydrodynamic Analysis of Marine Structures) is an open-source software for the numerical computation of the wave effect upon marine structures. It is based on boundary integral equations in the potential flow theory for analysis of wave-structure interactions. It is currently written in FORTRAN 90. The code has been developed by the author Yingyi Liu for almost a decade. 

HAMS is released in the hope that it will contribute to eliminating the inequality (for those who are not able to afford to purchase a costly commercial BEM software) in the continuous research developments related to offshore engineering and the ocean renewable energies.

HAMS is freely distributed under the Apache License, Version 2.0, http://www.apache.org/licenses/LICENSE-2.0, and may be modified and extended by researchers who intend to enhance its capabilities and port the code to other platforms. 

The success of HAMS should to a large portion be attributed to Prof. Bin Teng (Dalian University of Technology), who has tutored me the theory of potential flow in marine hydrodynamics and the programming skills using the [Boundary Element Method](https://en.wikipedia.org/wiki/Boundary_element_method). The code structure and the coding style of HAMS are exactly two of the examples that I have learned and inherited from Prof. Bin Teng.

## Theoretical Basis

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/global_coordinate_system.png" width="70%" /></p>
  
### - Please refer to the following papers for the theory:

The theory of panel method that has been used by HAMS is written in detail in the following two papers:

* Yingyi Liu (2019). "HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application." _Journal of Marine Science and Engineering_, 7: 81.

* Yingyi Liu, Changhong Hu, Makoto Sueyoshi, Hidetsugu Iwashita, Masashi Kashiwagi (2016). "Motion response prediction by hybrid panel-stick models for a semi-submersible with bracings." _Journal of Marine Science and Technology_, 21:742–757.

The deepwater Green function is using a fortran subroutine (https://github.com/Hui-Liang/Green-function-in-deep-water) developed by Dr. Hui Liang. For the detailed theory you may refer to the following three papers:

* Hui Liang, Huiyu Wu, and Francis Noblesse (2018). "Validation of a global approximation for wave diffraction-radiation in deep water." _Applied Ocean Research_, 74 : 80-86.

* Huiyu Wu, Hui Liang, and Francis Noblesse (2018). "Wave component in the Green function for diffraction radiation of regular water waves." _Applied Ocean Research_, 81: 72-75.

* Huiyu Wu, Chenliang Zhang, Yi Zhu, Wei Li, Decheng Wan, Francis Noblesse (2017). "A global approximation to the Green function for diffraction radiation of water waves." _European Journal of Mechanics-B/Fluids_, 65: 54-64.

The finite-depth Green function is using a fortran subroutine FinGreen3D (https://github.com/YingyiLiu/FinGreen3D) developed by Dr. Yingyi Liu. For the detailed theory you may refer to the following two papers:

* Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He (2018). "A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions." _Energy Conversion and Management_, 174: 516-536.

* Yingyi Liu, Hidetsugu Iwashita, Changhong Hu (2015). "A calculation method for finite depth free-surface green function." _International Journal of Naval Architecture and Ocean Engineering_, 7(2): 375-389.

Please cite appropriately the above papers in your relevant publications, reports, etc., if the HAMS code or its executable program has contributed to your work.

## Generated numerical results

### - Hydrodynamic coefficients

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/hydrodynamic_coefficients.png" width="60%" /></p>

### - Wave excitation force

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/wave_excitation_force_plot.png" width="65%" /></p>

### - Motion RAOs

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/RAO_of_motion.png" width="75%" /></p>

### - Free-surface elevation

  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/free_surface_elevation.png" width="60%" /></p>

## Features

### - Mesh element type

* HAMS can import meshes containing triangular panel type, quadrilateral panel type, or both.
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/panel_local_coordinates.png" width="90%" /></p>

### - OpenMP parallel processing

* HAMS can be run in the parallel mode using OpenMP techniques on PC's with multiple processors (CPU's).
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/OpenMP_parallel_process.png" width="95%" /></p>
  
### - Computational efficiency

* The following graph shows an example of DeepCwind semisubmersible using 8 threads for the computation:
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/computational_efficiency.png" width="80%" /></p>

## Useful Links

The following open-source software can be used to view the HAMS results: </p>
[1]. [BEMRosetta](https://github.com/izabala123/BEMRosetta). Developed by Iñaki Zabala, Markel Peñalba, Yerai Peña-Sanchez.<br/> 
[2]. [BEMIO](https://wec-sim.github.io/bemio/). Developed by National Renewable Energy Laboratory and Sandia 
Sandia National Laboratories. <br/> 

You may need HAMS to do the frequency-domain pre-processing before you use the following programs: </p>
[1]. [FAST](https://www.nrel.gov/wind/nwtc/fast.html) or [OpenFAST](https://openfast.readthedocs.io/en/master/). Developed by National Renewable Energy Laboratory.<br/> 
[2]. [WEC-Sim](https://wec-sim.github.io/WEC-Sim/). Developed by National Renewable Energy Laboratory and Sandia 
Sandia National Laboratories. <br/> 

## License

Code original author: Yingyi Liu (劉盈溢). [Google Scholar] (https://scholar.google.co.jp/citations?hl=ja&user=mpR3MvAAAAAJ&view_op=list_works&sortby=pubdate)

HAMS is free software: you can redistribute it and/or modify it under the terms of the Apache License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.

HAMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache License for details. You should have received a copy of the Apache License along with HAMS. If not, see http://www.apache.org/licenses/LICENSE-2.0 <br/>


