# OrbitalElements.jl

*Conversions between orbital constants for galactic dynamics in Julia.*

---
## Purpose and Scope

`OrbitalElements` provides changes of variables between various orbital constant/integrals in self-gravitating systems.
It is tailored for a specific set of systems, focusing on scenarios for linear response computations. 
By narrowing its scope, the library maintains clarity and efficiency in handling orbital dynamics within these systems. 
While it may not offer the generality of other related libraries ([galpy](https://www.galpy.org), [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama/tree/master)...), it excels in its targeted application domain.

For now, `OrbitalElements` specializes in handling orbits confined to a plane, such as in razor-thin discs or spherical systems, within central potentials. 
While its primary focus lies in these specific systems, it will be extended to encompass other systems in the future.

---
## Planar orbits

### Constants of motion

Given a central potential model, `OrbitalElements` provides mapping functions and jacobians 
between the following constants of motion
+ Model-independent constants
    - semi-major axis and eccentricity: $(a,e)$
    - pericentre and apocentre: $(r_p,r_a)$
+ "Analytic" constants
    - energy and angular momentum: $(E,L)$
+ "Integrated" constants
    - actions: $(J,L)$
    - frequency ratios: $(\alpha,\beta)$
    - frequencies: $(\Omega_1,\Omega_2)$

While [explicit formulae](formulae.md) for actions are often given as a function of $(r_p,r_a)$, the $(a,e)$-domain is more suitable to handle close-to-the-edge cases (near circular or near radial orbits) and to cure divergences in frequency computations.

Semi-major axis and eccentricity are therefore used as base constants, so that, e.g., going from actions to frequencies would in practice require performing $(J,L)\!\to\!(a,e)\!\to\!(\Omega_1,\Omega_2)$ (as no explicit formulae exists to generically convert actions in frequencies).



### Resonance variables

For a given potential, one can also compute the resonance mappings, called $(u,v)$, defined in [Fouvry & Prunet (2022)](https://doi.org/10.1093/mnras/stab3020). 
These mappings depend on the resonance vector $(n_1,n_2)$ and prove particularly useful for linear response' computations.
See `src/mappings/Resonance` for associated conversions.