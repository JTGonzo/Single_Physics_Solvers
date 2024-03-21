![](readme_support/xiao_white.gif#gh-light-mode-only)
![](readme_support/xiao_black.gif#gh-dark-mode-only)

# DISCLAIMER!!
These codes were very recently uploaded (March 20th, 2024) from my system and have yet to be vetted for general distribution/use. I am currently working through commenting, cleaning and re-testing these solver which were put together approximately 2 years ago. Please allow me the next few weeks to clean up the codes so as to make them more instructive, easy to use, and bug-free. 

`While this disclaimer is being displayed` **I strongly suggest the user refrain from using the 3D compressible flow solver** and exercise *some* caution when using the deconstructed **redbkit** solvers provided.

# CSM and CFD solver detailes
This repository contains two 2D and 3D high-fidelity computational fluid dyynamics (CFD) and computational structural mechanics solvers built (mainly) in MATLAB. Each solver irrespective of the domain is developed using a Galerkin least-squares finite element (FE) strategy. 

<ins>***3D Compressible Navier-Stokes Solver***</ins><br>
This positivity preserving variational Galerkin finite element scheme was developed at the [University of British Columbia](https://cml.mech.ubc.ca/) by the graduate students of my superviosor [Dr. Rajeev Jaiman](https://scholar.google.com/citations?user=iofAU68AAAAJ&hl=en&oi=ao). In this finite element solver an incompressible 3D Navier-Stokes solver is coupled in a partitioned manner to a homogenous mixture-based finite mass tranfer model that captures the cavitation physics. In addition, this cavitating flow solver is next coupled to a hybrid URANS-LES turbulence mode and temporal discretization is performed via a fully-implicit generalized-&alpha; time integration scheme. 

Further reading on this technique can be found in 
>[**[RKJ21] S. Kashyap and R. Jaiman. A robust and accurate finite element framework for cavitating flows with moving fluid-structure interfaces **, Computers and Mathematics with Applications, 2021.](https://doi.org/10.1016/j.camwa.2021.10.024)


<ins>***redbKIT Solver***</ins><br>
The redbKIT-based CFD and CSM solvers given are bare-bones deconstructions (/simplifications) of the original [toolkit](https://github.com/redbKIT/redbKIT). At face value they are exclusively the intellectual property of the original repository with a few original lines of functional code included for initialization, data output, and operational purposes. 

In the above solvers specifically (as opposed to the original toolkit), we exclusively use higher-order P<sub>2</sub> finite elements for the structural equations and stable Taylor-Hood P<sub>2</sub>/P<sub>1</sub> mixed finite elements for the fluid equations. A second-order backward difference scheme is adopted for temporal discretization in the CFD solvers whilst the generalized-&alpha; method is used for the CSM solvers. For the structural mechanics solvers (2D & 3D) the user is free to choose any one of three constitutive material models. They are; 
1. Linear Elastic 
2. St. Venant Kirchhoff 
3. Nearly-incompressible Neo-Hookean

As described in detail below; while the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit should be the primary reference point for the interested reader, it is my belief that the condensed solvers (as presented here in their simplified intermediary forms) hold their own unique value to the wider student/research community. In particular, for the purposes of streamlined learning, efficient tool adoption, and additional/alternative technique development. This was most certainly my case where the ability to quickly validate (/experiment with) the domain solvers themselves - when treated independently - made integrating them into a [partitioned fluid-structure interaction framework](https://github.com/JTGonzo/Partitioned_FSI) significantly easier. 

## Download, Installation and Running
-------

For detailed configuartion, installation and operation instructions, in particular for the redbKIT based solvers, please follow the steps contained in the brief [INSTALL.md](INSTALL.md) file provided.

## Further Clarifying Details
For absolute clarity, the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit by itself is an ***AMAZING*** resource for (young) academics seeking an education on the finite element method and reduced-order modelling. As such I **strongly** recommend the reader fork that repository and work through redbKIT's numerous well crafted example problems. *I personally give this resource the lion's share of credit in helping me complete my graduate studies.*

Having said that, however, as is true for (almost) all **well developed / robust** software packages used to tackle complex problems in computational physics; its sophistication tends to be its own undoing for those novice students beginning their journey in the simulation sciences or those young numerical methods researchers who are focusing on a niche element of the problem/framework/technique and just want the machinery around the issue to be (first) working and (secondly) easily digestible for peace of mind, debugging, development etc. purposes. As such, the codes provided above seek to shallow the learning curve and expedite the toolkit's adoption/integration. 

In the case of my graduate studies, for example, having simplified the redbKIT solvers into just those most essential numerical elements needed to perform reasonably complex CFD and CSM simulations, I was able to more seamless add (/test/troubleshoot) my own features to the already established physical solvers. In particular, the integration of conservative Dirchlet and Neummann interfacial boundary conditions into the fluid and structural solvers, respectively. These additions, as well as all associated interfacial data mapping and stabalization techniques also necessary, then allowed me to developm of my own [partitioned numerical framework](https://github.com/JTGonzo/Partitioned_FSI) which I used to investigate ["*coupling  instability in low mass-ratio partitioned FSI simulations*"](https://jtgonzo.github.io/).  

## redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a MATLAB library for finite element simulation and reduced-order modeling of Partial Differential Equations developed at [EPFL](https://www.epfl.ch/) - [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)). 

In particular, it includes straightforward implementations of many of the algorithms presented in the companion book:

>[**[QMN16] A. Quarteroni, A. Manzoni, F. Negri. Reduced Basis Methods for Partial Differential Equations. An Introduction**, Springer, 2016.](http://www.springer.com/us/book/9783319154305#aboutBook)

For further details please visit [redbKIT website](http://redbkit.github.io/redbKIT/).

### License
-------

**redbKIT is distributed under BSD 2-clause license**

Copyright (c) 2015-2017, Ecole Polytechnique Fédérale de Lausanne (EPFL)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### **redbkit** Development
-------

redbKIT was developed and maintained by [`Federico Negri`](https://www.linkedin.com/in/negrifederico/). Paola Gervasio (Università degli Studi di Brescia) is gratefully acknowledged for granting the use of parts of the finite element code MLife.


## Contact
-------
If you have any questions regarding the contents of this repository or the use of these solvers, please feel to email me at <josetg@vt.edu>.
