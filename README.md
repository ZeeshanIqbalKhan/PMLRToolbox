# PMLR Toolbox
Official repository of paper, entitled ***"Nonlinear Control Allocation Using A Piecewise Multi-Linear Representation"***.

[https://arxiv.org/abs/2208.10411](https://arxiv.org/abs/2208.10411)

## Abstract
Nonlinear control allocation is an important part of modern nonlinear dynamic inversion based flight control systems which require highly accurate model of aircraft aerodynamics. Generally, an accurately implemented onboard model determines how well the system nonlinearities can be canceled. Thus, more accurate model results in better cancellation, leading to the higher performance of the controller. In this paper, a new control system is presented that combines nonlinear dynamic inversion with a piecewise multi-linear representation based control allocation. The piecewise multi-linear representation is developed through a new generalization of Kronecker product for block matrices, combined with the canonical piecewise linear representation of nonlinear functions. Analytical expressions for the Jacobian of the piecewise multi-linear model are also presented. Proposed formulation gives an exact representation of piecewise multi-linear aerodynamic data and thus is capable of accurately modeling nonlinear aerodynamics over the entire flight envelope of an aircraft. Resulting nonlinear controller is applied to control of a tailless flying wing aircraft with ten independently operating control surfaces. The simulation results for two innovative control surface configurations indicate that perfect control allocation performance can be achieved, leading to better tracking performance compared with ordinary polynomial-based control allocation. 


If you use this work for your projects, please take the time to cite our paper:

  ```
  @misc{rajput2022,
      title = {Nonlinear Control Allocation Using A Piecewise Multi-Linear Representation},
      author = {Rajput, Jahanzeb and Khan, Hafiz Zeeshan Iqbal},
      year = {2022},
      eprint={2208.10411},
      archivePrefix={arXiv},
      primaryClass={eess.SY}
  }
  ```
