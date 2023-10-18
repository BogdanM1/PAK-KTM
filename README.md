# PAK-KTM
Description:
PAK-KTM is high-performance software for finite element analysis (FEA), developed at the University of Kragujevac and the Research and Development Center for Bioengineering, BioIRC. 
The program is written in FORTRAN 77/90/95 and is capable of solving diffusion and convection fields, enabling the modeling of mass transport of ions or molecules.
The fundamental laws utilized include Fick's law of diffusion, from which the mass balance equation is derived and translated into a form applicable in the finite element method. 
The model incorporates the assumption that the transport of molecules can be described by a nonlinear, non-stationary diffusion process.

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


Below is the description in Serbian.

Опис:
PAK-KTM је софтвер високих перформанси за анализу методом коначних елемената (МКЕ), развијен на Универзитету у Крагујевцу и у Истраживачко-развојном центру за биоинжињеринг, БиоИРЦ, 
Програм је развијен у FORTRAN-у 77/90/95. PAK-KTM решава дифузионо и конвекционо поље, тј. омогућава моделирање транпорта масе јона или молекула. 
Основни закони су Фиков закон дифузије, за кога је даље изведена једначина баланса масе и преведена у облик применљив у методи коначних елемената.
У моделу је уведена претпоставка да се транспорт молекула може описати нелинеарно-нестационарним дифузионим процесом. 
Други закон је Дарсијев закон који описују транспорт флуида у порозним срединама применом конвекције.

<img src="imgs/Diag1.png?raw=true" title="Ток података" width="200" align="center">

<img src="imgs/Diag2.png?raw=true" title="Којић транспорт елемент" width="100" align="center">

<img src="imgs/Diag3.png?raw=true title="Конективни елементи" width="100" align="center">

<img src="imgs/results1.png?raw=true title="Поређење smeared и детаљног модела" width="100" align="center">






