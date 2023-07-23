# PB-LJ EDL Model

In this code, we demonstrate the solution of the theory developed in <a href="http://dx.doi.org/10.1039/D3CP00745F">Incorporating Ion-Specific van der Waals and Soft Repulsive Interactions in the Poisson-Boltzmann Theory of Electrical Double Layers</a> to incorporate van der Waals and soft repulsive interactions into a Poisson-Boltzmann framework. There are five MATLAB codes, as documented below:

1. ion_wall_code: We solve for the potential and ionic concentration profiles inside an electrical double layer (EDL) with ion-wall interactions included in the model. Note that if you want to use the ion-wall PB-LJ code without the approximation ($u^{iw}_{LJ}(L/2) ≈ 0$), use the version in the Comparison_Aluru folder.

2. ion_ion_code_both.m: In this code, we solve for the potential and ionic concentration profiles inside an EDL with ion-ion and ion-wall interactions included in the model.

3. Capacitance_NaF_Ag_111.m: In this code, we plot differential capacitance curves for aqueous NaF solutions in contact with an Ag(111) electrode and compare the results with experimental data (Valette, J. Electroanal. Chem. Interfacial Electrochem.,1989, 269, 191–203).

4. Comparison_Aluru: In this code, we compute concentration profiles using our PB-LJ model to compare with molecular dynamics results (Mashayak and Aluru, J. Chem. Phys., 2017, 146, 044108).

6. sigma_epsilon_fitting.m: In this code, we fit the Lennard-Jones (LJ) parameters for the polyatomic sulfate ion using the LJ parameters for the constituent sulfur and oxygen atoms.

## Citation

```bibtex
@article{D3CP00745F,
title  ="Incorporating Ion-Specific van der Waals and Soft Repulsive Interactions in the Poisson-Boltzmann Theory of Electrical Double Layers",
author ="Seal, Aniruddha and Tiwari, Utkarsh and Gupta, Ankur and Govind Rajan, Ananth",
journal  ="Phys. Chem. Chem. Phys.",
year  ="2023",
pages  ="-",
publisher  ="The Royal Society of Chemistry"
}
```
