<h1 align="center"> PAK-KTM</h1>

<p align="justify">
 PAK-KTM (Program for Structural Analysis - *Kojic Transport Model*) is high-performance software for finite element analysis (FEA), developed at the University of Kragujevac and the Research and Development Center for Bioengineering, BioIRC. 
The program is written in FORTRAN 77/90/95 and is capable of solving diffusion and convection fields, enabling the modeling of mass transport of ions or molecules.
The fundamental laws utilized include Fick's law of diffusion, from which the mass balance equation is derived and translated into a form applicable in the finite element method. 
The model incorporates the assumption that the transport of molecules can be described by a nonlinear, non-stationary diffusion process.

 The software involves modeling systems for the analysis and assessment of convective-diffusive transport in complex systems such as tumors and organs consisting of a large number of capillaries and subdomains (capillary domain, lymphatic domain, tissue domain, cell domain, intercellular space domain, etc.). Determination of basic mass transport process variables such as velocities, pressure, or concentration is performed using the finite element method (FEM). 

Modern medicine involves the application of the latest technologies for more efficient human treatment and predicting disease development, as well as predicting drug delivery within complex biological systems. As a result of numerical calculations using this method, significant information is obtained, along with a better understanding of transport processes occurring in the vascular system, organs, and intercellular space. Users of this technical solution can model complex tissue and organ systems in a very simple way, using distributed (**smeared**) fields, significantly facilitating the generation of complex systems and organs and further reducing costs for pre-clinical and clinical trials.
</p>

Another law employed is Darcy's law, which describes the transport of fluids in porous media through convection.
When using this software, please cite following papers:
- M. Kojic, M. Milosevic, N. Kojic, Z. Starosolski, K. Ghaghada, R. Serda, A. Annapragada, M. Ferrari, and A. Ziemys, “A multi-scale FE model for convective–diffusive drug transport within tumor and large vascular networks,” Computer Methods in Applied Mechanics and Engineering, vol. 294. Elsevier BV, pp. 100–122, Sep-2015. [http://dx.doi.org/10.1016/j.cma.2015.06.002]

- M. Kojic, M. Milosevic, N. Kojic, E. J. Koay, J. B. Fleming, M. Ferrari, and A. Ziemys, “Mass release curves as the constitutive curves for modeling diffusive transport within biological tissue,” Computers in Biology and Medicine, vol. 92. Elsevier BV, pp. 156–167, Jan-2018. [https://doi.org/10.1016/j.compbiomed.2016.06.026]

- M. Kojic, M. Milosevic, V. Simic, E. J. Koay, J. B. Fleming, S. Nizzero, N. Kojic, A. Ziemys, and M. Ferrari, “A composite smeared finite element for mass transport in capillary systems and biological tissue,” Computer Methods in Applied Mechanics and Engineering, vol. 324. Elsevier BV, pp. 413–437, Sep-2017. [http://dx.doi.org/10.1016/j.cma.2017.06.019]

- M. Kojić, M. Milošević, V. Simić, E. J. Koay, N. Kojić, A. Ziemys, and M. Ferrari, “EXTENSION OF THE COMPOSITE SMEARED FINITE ELEMENT (CSFE) TO INCLUDE LYMPHATIC SYSTEM IN MODELING MASS TRANSPORT IN CAPILLARY SYSTEMS AND BIOLOGICAL TISSUE,” Journal of the Serbian Society for Computational Mechanics, vol. 11, no. 2. Faculty of Engineering, University of Kragujevac, pp. 108–119, Dec-2017. [http://doi.org/10.24874/jsscm.2017.11.02.09]

- M. Milosevic, V. Simic, B. Milicevic, E. J. Koay, M. Ferrari, A. Ziemys, and M. Kojic, “Correction function for accuracy improvement of the Composite Smeared Finite Element for diffusive transport in biological tissue systems,” Computer Methods in Applied Mechanics and Engineering, vol. 338. Elsevier BV, pp. 97–116, Aug-2018.[https://doi.org/10.1016/j.cma.2018.04.012]

- M. Kojic, M. Milosevic, V. Simic, E. J. Koay, N. Kojic, A. Ziemys, and M. Ferrari, “Multiscale smeared finite element model for mass transport in biological tissue: From blood vessels to cells and cellular organelles,” Computers in Biology and Medicine, vol. 99. Elsevier BV, pp. 7–23, Aug-2018. [https://doi.org/10.1016/j.compbiomed.2018.05.022]

- M. Milosevic, D. Stojanovic, V. Simic, B. Milicevic, A. Radisavljevic, P. Uskokovic, and M. Kojic, “A Computational Model for Drug Release from PLGA Implant,” Materials, vol. 11, no. 12. MDPI AG, p. 2416, 29-Nov-2018. [http://dx.doi.org/10.3390/ma11122416]

- R. Santagiuliana, M. Milosevic, B. Milicevic, G. Sciumè, V. Simic, A. Ziemys, M. Kojic, and B. A. Schrefler, “Coupling tumor growth and bio distribution models,” Biomedical Microdevices, vol. 21, no. 2. Springer Science and Business Media LLC, 25-Mar-2019. [https://doi.org/10.1007/s10544-019-0368-y]


<p align="center">
 <figure>
<img src="imgs/Diag1.png?raw=true"  width="300" />
  <figcaption>Data flow</figcaption>
 </figure>

 <figure>
<img src="imgs/Diag22.PNG?raw=true"  width="300" />
  <figcaption>Kojic transport element</figcaption>
 </figure>

 <figure>
<img src="imgs/Diag3.png?raw=true" width="300" />
  <figcaption>Connective elements</figcaption>
 </figure>

</p>


<h2 align="center"> Requirements</h2>

* **Windows:** <br/>
* VisualStudio >=2017
* Intel® oneAPI Base Toolkit
* Intel® Distribution for GDB*
* Intel® oneAPI HPC Toolkit
  <br/>
  <br/>
  
* **Linux:**<br/>
* GNU Fortran (GCC) 4.4.7 20120313
* The C compiler GNU 4.4.7
* MUMPS 5.0.2.

<h2 align="center"> Manual (Windows) </h2>
Download code from github and open visual studio. As start-up project set PAKKTM and build project.
Once the project is succesfully built you can run pakktm by pressing the Start button as shown below.

<p align="center">
<img src="imgs/Manual1.PNG?raw=true" width="300" />
</p>

<p align="center">
<img src="imgs/Manual2.PNG?raw=true" width="300" />
</p>

<h2 align="center"> Manual (Linux) </h2>
Download code from github navigate to pakktm/build and run the following commands:<br />
<p align="left">
cmake .. <br />
make <br />
./pakktm <br />
</p>
The program will prompt you to type-in the input dat file. <br />

<h2 align="center"> Examples </h2>

<h3 align="center"> Mouse liver model </h3>

<p align="justify">
The finite element model consists of 1D pipe FEs for larger vessels (7736 elements), 3D composite smeared elements (39832 elements), and connectivity elements (726 elements) for connecting large vessels with continuum nodes (capillary domain DOF) of smeared FEs. There are two separate tumors within liver, with a total of 316 elements. The total number of nodes is 54590.
</p>

<p align="center">
 <figure>
<img src="imgs/examples_imgs/liver/1.jpg?raw=true"  width="300" />
  <figcaption>Liver model: geometry, tumor domains and pressures within large vessels</figcaption>
 </figure>

 <figure>
<img src="imgs/examples_imgs/liver/2.jpg?raw=true"  width="300" />
  <figcaption>Concentration field in liver with tumors (marked with dashed lines), dotted results in tissue domain and with full mesh in tumors, for times t =10, 20 and 50s.</figcaption>
 </figure>

 <figure>
<img src="imgs/examples_imgs/liver/3.jpg?raw=true" width="300" />
  <figcaption>Pressure fields for two views: a) Full mesh; b) Clipped mesh; c) Dotted representation of results in tissue, and with full mesh in tumors.</figcaption>
 </figure>

</p>


<h3 align="center"> Pancreas model </h3>


<h3 align="center"> PLGA </h3>









