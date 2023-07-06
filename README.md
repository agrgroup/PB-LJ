# PB-LJ EDL Model

In this code, we demonstrate the solution of the theory developed in <a href="https://arxiv.org/abs/2302.07628">Incorporating Ion-Specific van der Waals and Soft Repulsive Interactions in the Poisson-Boltzmann Theory of Electrical Double Layers</a> to incorporate van der Waals and soft repulsive interactions into a Poisson-Boltzmann framework. There are four MATLAB codes, as documented below:

1. ion_wall_code: We solve for the potential and ionic concentration profiles inside an electrical double layer (EDL) with ion-wall interactions included in the model.

2. ion_ion_code_both.m: In this code, we solve for the potential and ionic concentration profiles inside an EDL with ion-ion and ion-wall interactions included in the model.

3. Capacitance_NaF_Ag_111.m: In this code, we plot differential capacitance curves for aqueous NaF solutions in contct with a Ag(111) electrode, and compare the results with experimental data.

4. sigma_epsilon_fitting.m: In this code, we fit the Lennard-Jones (LJ) parameters for the polyatomic sulfate ion using the LJ parameters for the constituent sulfur and oxygen atoms.

## Citation

```bibtex
@article{seal2023incorporating,
  title={Incorporating Ion-Specific van der Waals and Soft Repulsive Interactions in the Poisson-Boltzmann Theory of Electrical Double Layers},
  author={Seal, Aniruddha and Tiwari, Utkarsh and Gupta, Ankur and Rajan, Ananth Govind},
  journal={arXiv preprint arXiv:2302.07628},
  year={2023}
}
```
