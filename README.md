# OpenPhase4PVD
Multi-phase field model to simulate the PVD, achieved by OpenPhase

## modelGen1
The model is based on the multi-phase field model developed by Ingo Steinbach et al. and OpenPhase(https://github.com/ICAMS/OpenPhase).
Here, the deposition flux is coupled in the model by merging a driving force term of deposition. So firstly the kinetic deposition flux ought to be translated to the thermodynamic driving force term.
Phase fields are only considered in this model.

The method to run the simulation:
''' python
python simulate.py > simulation.log
'''

## modelGen2
Based on modelGen1, a mass field $\rho$ is introduced in modelGen2.
And the governing equation for mass field would be solved:
$$
d rho/dt = div(D_phi_v * grad(rho)) - div(phi_v * rho * v)) - C_dep * max(d phi_s/dt, 0)
$$
A vapor transport is deposited from the mother phase to the solid surface.
