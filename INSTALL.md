===========================================
Instructions - rebkit-based solvers 
===========================================

Pre-requisites
--------------

For the basic installation, you need:<br>
-> a reasonably recent version of Matlab (anything above R2017a should be fine)<br>
-> a supported C compiler, see e.g. http://it.mathworks.com/support/compilers/R2016a/<br>
   (gcc on Unix and gcc/clang on Mac OS X work smoothly)


Basic Installation
------------------

Open Matlab and navigate to the root folder of the downloaded solver that you will be using (e.g. **3D_structural_mechanics**). Then type 

>> make

to compile the C-assembly routines and "mexify" some other files.


Running the Example Provided 
------------------

With Matlab opened; in the same folder as the *make.m* file ran previously, type 

>> main_CFD

OR 

>> main_CSM

as determined by the solver type. Doing this begins the example simulation and creates a "Figures" folder where the simulation configuration is output periodically as a series of *.vtk files. In the casse of the CFD solver a "Results" folder is also created where the resultant aerodynamic forces (integrated over the defined FSI interface) are print to a file contained therein after each time step. 

Running your own Simulation
------------------
<ins>***Mesh and .MAT file generation***</ins><br>
See this brief tutorial on how to set up the **.mat** files needed to perform your own simulations. These .mat files define the geometric/spatial properties of the problem you intend to investigate. Pay special attention to how you define and identify the boundary and interface surfaces. 

<ins>***Editting main_CFD.m or main_CSM.m***</ins><br>
The only adjustments that are needed in this file to run your own simulations are 

*line 8 :* load('beam.mat')

*line 15 :* vtk_filename = 'Figures/Flap_';

These lines are where you indicate the names you chose for the domain specific .mat files that you will be importing into the solver as well as the prefix string that you assign to your output data.  

<ins>***Editting CFD_data.m or CSM_data.m***</ins><br>
These two files are the primary locations where you define all the domain specific parameters needed for you unqiue problem. Parameters that you are able to set include: 

1. material properties 
2. boundary conditions
3. initial conditions
4. structure's constitutive model to use
5. temporal and spatially varying properties
6. time step/integration properties
7. aerodynamic forces output variables


===========================================
Instructions - multi=phase CFD solver 
===========================================

Pre-requisites
--------------

For the basic installation, you just need a reasonably recent version of Matlab (anything above R2015a should be fine). There are no additional compiler or system requirements. 


Running the Example Provided 
------------------

Open Matlab and navigate to the root folder of the downloaded solver that you will be using (e.g. **3D_compressible_NS**). Then type 

>> main_CFD

doing this begins the example simulation provided. 

