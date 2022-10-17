
# OrbitalElements.jl

`OrbitalElements.jl` is a package written in Julia to compute numerical elements for astronomical orbits to high precision, for arbitrary potentials.

-----------------------------

### Quick activate

`OrbitalElements` is (currently) unregistered, and as such if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Shortcut version: after cloning the repository, navigate to the top directory. Start julia (`julia`), then the package manager (`]`), then register the package in development mode (`dev .`).

Then, in your program, if you want to access specific elements listed below, use `import OrbitalElements`.

-----------------------------

### Sample potentials

`OrbitalElements` ships with a handful of simple potentials and distribution functions for testing purposes. These are located in `src/Potential`: (1) the 3d Isochrone potential, (2) the 3d Plummer potential, (3) the 2d Mestel-Zang disc, and (4) the 2d Kuzmin-Toomre disc. Corresponding distribution functions are located in `src/DistributionFunctions`.

----------------------------

### Obtaining Orbital Frequencies

The principal use for `OrbitalElements` is to provide descriptions of orbits. This means
1. Conversion from (semimajor axis, eccentricity) to (pericentre, apocentre) to (energy, angular momentum), to coordinates aligned with resonance vectors.
2. Additional support for actions and angles.

`ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)` will compute frequencies (Ω₁,Ω₂) given a potential (ψ) plus two derivatives (dψ,d2ψ), for an orbit described by semimajor axis and eccentricity.

`Θ(ψ,dψ,d2ψ,u,rp,ra)` will compute the anomaly along an orbit, given a potential plus two derivatives.

-----------------------------

### Mapping to Resonance Space

For a given potential, one can also compute the resonance mappings, called (u,v). See `src/Resonance` for associated conversions.
`UVFromαβ(α,β,n₁,n₂,ωmin,ωmax)` will compute the resonant mappings (u,v) for a given frequency pair (α,β), resonance vector (n₁,n₂), and frequency limits.

-----------------------------

### Notes
By default, `OrbitalElements` uses pericentre and apocentre. If you want to use semimajor axis and eccentricity units, transformations are available. `RpRaFromAE(a,e)` will return pericentre and apocentre from semimajor axis and eccentricity.

`OrbitalElements` also uses the Henon (1971) technique to cure radial velocity divergences at peri- and apocentre. One could use other methods; the software is constructed to allow for drop-in anomaly replacements. See `src/Henon`.

-----------------------------

### Authors

Mike Petersen -  @michael-petersen - petersen@iap.fr, michael.petersen@roe.ac.uk

Mathieu Roule -  @MathieuRoule     - roule@iap.fr
