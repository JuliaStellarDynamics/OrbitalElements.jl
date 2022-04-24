
# OrbitalElements.jl

`OrbitalElements.jl` is a package written in Julia to compute numerical elements for astronomical orbits..

-----------------------------

## Quick activate

As `OrbitalElements` is unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/JuliaOrbitElements.git"`. If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.

Another option, for development: after cloning the repository, in the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using OrbitalElements` into the Julia interpreter. You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import OrbitalElements` which will give access modeled on `OrbitalElements.rpra_from_ae` (for example).

If you need to remove the package for some reason (force installing a new version, for instance), you can `import Pkg;Pkg.rm("OrbitalElements")`.


-----------------------------

### Obtaining Orbital Frequencies

`compute_frequencies(ψ,dψ/dr,d²ψ/dr²,rp,ra)` will compute frequencies (Ω₁,Ω₂) given a potential (ψ) plus two derivatives (dψ/dr,d²ψ/dr²), for an orbit described by pericentre (rp) and apocentre (ra).

Example for matching against isochrone:
```
using OrbitalElements
rp,ra = 0.9,1.1
OrbitalElements.compute_frequencies_rpra(OrbitalElements.isochrone_psi,OrbitalElements.isochrone_dpsi_dr,OrbitalElements.isochrone_ddpsi_ddr,rp,ra)
```

`Theta(ψ,dψ/dr,d²ψ/dr²,u,rp,ra)` will compute the anomaly along an orbit, given a potential plus two derivatives.

-----------------------------

### Mapping to Resonance Space

For a given potential, one can also compute the resonance mappings, called (u,v).
`uv_from_alphabeta(α,β,n₁,n₂,dψ/dr,d²ψ/dr²)` will compute the resonant mappings (u,v) for a given frequency pair (α,β), resonance vector (n₁,n₂), and potential derivatives.

-----------------------------

### Analytic isochrone quantities for matching
We know analytic isochrone frequencies and actions. The definitions are included as a verification of the methods.

`isochrone_Omega_1_2(rp,ra,bc,M,G)` will return the analytic isochrone frequencies (Ω₁,Ω₂) for an orbit defined by pericentre (rp), apocentre (ra) and potential sourced by (bc,M,G; all unity by default).

-----------------------------

### Notes
By default, `OrbitalElements` uses pericentre and apocentre. If you want to use semimajor axis and eccentricity units, transformations are available. `rpra_from_ae(a,e)` will return pericentre and apocentre from semimajor axis and eccentricity.

-----------------------------

### Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
