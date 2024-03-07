---
# Potentials

`OrbitalElements` ships with a handful of simple potential models for testing purposes. These are located in `src/Potential`: (1) the 3d Isochrone potential, (2) the 3d Plummer potential, (3) the 2d Mestel-Zang disc, and (4) the 2d Kuzmin-Toomre disc. Corresponding distribution functions are located in `src/DistributionFunctions`.

---
## Obtaining Orbital Frequencies

The principal use for `OrbitalElements` is to provide descriptions of orbits. This means
1. Conversion from (semimajor axis, eccentricity) to (pericentre, apocentre) to (energy, angular momentum), to coordinates aligned with resonance vectors.
2. Additional support for actions and angles.

`frequencies_from_ae(a,e,model,params)` will compute frequencies (Ω₁,Ω₂) given a potential model, for an orbit described by semimajor axis and eccentricity.

---
## Obtaining resonant coordinates
`UVFromαβ(α,β,n₁,n₂,ωmin,ωmax)` will compute the resonant mappings (u,v) for a given frequency pair (α,β), resonance vector (n₁,n₂), and frequency limits.