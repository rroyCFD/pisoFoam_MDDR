## Application : pisoMoukalledFoam (piso algorithm with SIMPLE & PRIME iterators)

### Authors

- Luca Mangani   HSLU email: luca.mangani@hslu.ch
- Marwan Darwish AUB  email: darwish@aub.edu.lb
- Fadl Moukalled AUB  email: moukalled@aub.edu.lb

#### Modification and adaptation for OF-6
Rajib Roy  University of Wyoming email: rroy@uwyo.edu, roy.rajib@live.com

### Description
Transient solver for incompressible, turbulent flow, using the PISO algorithm. The PISO algorithm followes the original paper by Issa:

    Issa RI (1982) Solution of the implicit discretized fluid flow equations by operator splitting. Mechanical Engineering Report, FS/82/15, Imperial College, London,

and later modification proposed by Yen et al.

    Yen RH, Liu CH (1993) Enhancement of the SIMPLE algorithm by an additional explicit corrector step. Numer Heat Transfer, Part B 24:127–141

The PISO algorithm is illustrated in the book [The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM® and Matlab](https://www.springer.com/gp/book/9783319168739) by  Moukalled, F., Mangani, L., Darwish, M. in section _15.7.3_.


    1. To compute the solution at time t + Δt, use as an initial guess the solution at time t for pressure, velocity, and mass flow rate fields p(n), u(n), and phi(n) respectively.

** SIMPLE Steps **

    2. Solve implicitly the momentum equation given by Eq. (15.70) to obtain a new velocity field v*.

    3. Update the mass flow rate at the cell faces using the Rhie-Chow interpolation technique (Eq. 15.100) to obtain a momentum satisfying mass flow rate field phi* .

    4. Using the new mass flow rates, assemble the pressure correction equation (Eq. 15.98) and solve it to obtain a pressure correction field p′.

    5. Update the pressure and velocity fields at the cell centroids and the mass flow rate at the cell faces to obtain continuity-satisfying fields using Eq. (15.101).

** PRIME Steps ** 

    6. Using the latest available velocity and pressure fields, calculate the coefficients of the momentum equation and solve it explicitly.
    
    7. Update the mass flow rate at the cell faces using the Rhie-Chow interpolation technique.

    8. Using the new mass flow rates, assemble the pressure correction equation (Eq. 15.183) and solve it to obtain a pressure correction field.

    9. Update pressure, velocity, and mass flow rate fields using expressions similar to the ones given in Eq. (15.101).

    10. Go to step 6 and repeat based on the desired number of corrector steps.
    
    11. Set the initial guess for velocity, mass flow rate, and pressure as u ** , phi**, and p* .

    12. Go back to step 2 and repeat until convergence.

    13. Set the solution at time t + Δt to be equal to the converged solution and set the current time t + Δt to be t.

    14. Advance to the next time step.

    15. Go back to step 1 and repeat until the last time step is reached.



### Disclaiimer:

This application is built based on [OpenFOAM version-6](https://openfoam.org/release/6/). Please read the _About OpenFOAM_ section to learn more on OpenFOAM.

The application is free to use. The author neither provide any warranty nor shall be liable for any damage incurred from this application.


#### About OpenFOAM

OpenFOAM is the leading free, open source software for computational fluid dynamics (CFD), owned by the OpenFOAM Foundation and distributed exclusively under the [General Public Licence (GPL)](http://www.gnu.org/copyleft/gpl.html). The GPL gives users the freedom to modify and redistribute the software and a guarantee of continued free use, within the terms of the licence. To learn more visit [https://openfoam.org/](https://openfoam.org/)
