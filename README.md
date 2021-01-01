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
### - please refer to the following papers:

Yingyi Liu (2019). HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application. Journal of Marine Science and Engineering 7: 81.

Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He (2018). A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions. Energy Conversion and Management, 174: 516-536.

Please cite the above papers in your relevant publications if the HAMS code or its executable program has contributed to your work.

## Features
### - Supported file formats

* Hydrodynamic coefficients
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/hydrodynamic_coefficients.png" width="1000" title="Hydrodynamic coefficients"></p>
* Wave excitation force
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/wave_excitation_force_plot.png" width="800" title="Wave excitation force"></p>
* Motion RAOs
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/RAO_of_motion.png" width="700" title="Motion RAOs"></p>
* Free-surface elevation
  <p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/free_surface_elevation.png" width="600" title="Free-surface elevation"></p>

### - Load the hydrodynamic coefficients from one format and save them in another

The goal is to have a good robustness in the handling of files


### - Compare the hydrodynamic coefficients for the same geometry from different software

- Damping for the same geometry got from different solvers
  
<p align="center"><img src="https://github.com/YingyiLiu/HAMS/blob/master/Other/md_resources/hydrodynamic_coefficients.png" width="300" title="Damping for the same geometry got from different solvers"></p>

- Excitation force for the same geometry got from different solvers_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/2%20solvers%20exc.jpg" width="300" title="Excitation force for the same geometry got from different solvers"></p>

### - Forces handling

It simmetrizes the available forces in all directions, averaging them when they are available on both possitive and negative headings. Some examples cases:
* Only the forces on positive headings from 0 to 180º have been processed: Symmetrize duplicates them to the negative heading values from 0 to -180º
* Both positive and negative headings forces have been processed: Symmetrize averages them

### - Compare the hydrodynamic coefficients for the same geometry for different discretization levels
### - Compare the hydrodynamic coefficients for different geometries

- Damping for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20B.png" width="300" title="Damping for different offshore wind floating platforms"></p>

- Excitation force for different offshore wind floating platforms_
  
<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/offshore%20wind%20platforms%20exc.jpg" width="300" title="Excitation force for different offshore wind floating platforms"></p>

### - FOAMM connection

[Finite Order Approximation by Moment-Matching (FOAMM)](http://www.eeng.nuim.ie/coer/wp-content/uploads/2019/02/FOAMM-Manual.pdf) is an application developed by N. Faedo, Y. Peña-Sanchez and J. V. Ringwood in the [Maynooth University](https://www.maynoothuniversity.ie/)'s [Centre for Ocean Energy Research (COER)](http://www.eeng.nuim.ie/coer/), that implements the moment-matching based frequency-domain identification algorithm.

BEMRosetta allows an interactive and seamless FOAMM connection to get state space coefficients.

### - Mesh loading, combining them for visual comparison 

Several meshes can be loaded in this basic viewer, allowing a visual comparison of geometries.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/deepcwind.jpg" width="300" title="Mesh loading"></p>


### - Mesh handling

- Interactive mesh rotation and translation around user defined center
- Automatic free surface, underwater surface, center of buoyancy, hydrostatic stiffness matrix, and other parameters calculation
- Improved viewer including dropdown menu in viewer screen
- Hydrostatic stiffness matrix viewer
- Mesh healing option
    
### - Nemoh

Added Nemoh launcher. It can load an existing Nemoh.cal file, lets you editing it, and creates the set of files to launch Nemoh from a .bat file (it replaces the classic MATLAB launcher)

### - Bonus

BEMRosetta now includes an OpenFAST .out/.outb file reader, designed to be very easy to use.
It includes an online mode for updating the plot while the simulation is in progress.
Files may be opened by drag and drop, and parameters are filtered by name or units.

<p align="center"><img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/FAST_Reader.png" width="800" title="OpenFAST .out/.outb reader"></p>

### - Other

All files, mesh, Nemoh or BEM files, can be loaded by Drag and Drop or Copy and Paste from file explorer in Windows and Linux.

<p align="center">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Underwater.png" width="45%" title="Underwater mesh and waterline">
  <img src="https://github.com/izabala123/BEMRosetta/blob/master/other/md%20resources/Mesh.png" width="45%" title="All mesh ans waterline">
</p>


## Useful Links

The following open-source software can be used to view the HAMS results: <br/>
[1].[BEMRosetta](https://github.com/izabala123/BEMRosetta).<br/>
[2].[BEMIO](https://wec-sim.github.io/bemio/).<br/>


## License

HAMS is free software: you can redistribute it and/or modify it under the terms of the Apache License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.\

HAMS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for details. You should have received a copy of the GNU General Public License along with BEMRosetta. If not, see http://www.gnu.org/licenses/.<br/>


