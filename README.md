
# OrbitalElements.jl

`OrbitalElements.jl` is a package written in Julia to compute numerical elements for astronomical orbits..

-----------------------------

## Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using OrbitalElements` into the Julia interpreter. You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import OrbitalElements` which will give access modeled on `OrbitalElements.rpra_from_ae` (for example).

As `OrbitalElements` is unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

-----------------------------

## Obtaining Orbital Frequencies

`compute_frequencies_henon(potential,dpotential,ddpotential,rp,ra)` will compute frequencies give a potential plus two derivatives, for an orbit described by pericentre (rp) and apocentre (ra)

Example for matching against isochrone:
`compute_frequencies_henon(OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr,0.9,1.1)`

-----------------------------

## Matching Isochrone Frequencies
We know analytic isochrone frequencies and actions, so this makes sense for a test case.

`isochrone_Omega_1_2(rp,ra,bc,M,G)` will return the analytic isochrone frequencies for an orbit defined by pericentre (rp), apocentre (ra) and potential sourced by (bc,M,G; all unity by default).

-----------------------------

## Notes
By default, `OrbitalElements` uses pericentre and apocentre. If you want to use semimajor axis and eccentricity units, transformations are available.

`rpra_from_ae(a,e)` will return pericentre and apocentre from semimajor axis and eccentricity.

-----------------------------

## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
