---
title: 'BlackHoleImages: A Wolfram Mathematica paclet for analytical black hole imaging'
tags:
  - Accretion Disks
  - Analytical Methods
  - Black Holes
  - Relativistic Astrophysics
  - Wolfram Mathematica
authors:
  - name: David Podrapsky
    orcid: 0009-0000-7558-2509
    affiliation: 1
  - name: Vojtech Witzany
    orcid: https://orcid.org/0000-0002-9209-5355
    affiliation: 1
affiliations:
  - name: Institute of Theoretical Physics, Charles University, Czechia
    index: 1
date: 16 February 2025
bibliography: paper.bib

---

# Summary

Black holes are astrophysical objects for which the gravitational field is strong enough that even particles moving at the speed of light cannot escape their influence. For their proper visualization, general relativistic effects on photon trajectories must thus be implemented. These include the changed shape of the trajectories, but also the changes in frequency and intensity, which is especially of interest when we study accretion disks forming in the black hole's proximity.

# Statement of need

`BlackHoleImages` is a Wolfram Mathematica paclet for generating images of accretion disks or stellar background in the Kerr geometry, which is usable for rotating, uncharged black holes. It contains ready-to-use functions in the `KerrImages` package, which allow imaging of equatorial disks described by the simple "alpha disk" model proposed by Shakura and Sunyaev [@Shakura], and also the stellar background distorted by the Kerr geometry.

The disk model is implemented in the `AlphaDiskModel` package, which can be redesigned as needed in future implementations. The analytical implementation of null geodesics in the Kerr spacetime is contained within the `KerrNullGeodesics` package. The solutions themselves were mostly taken from Gralla and Lupsasca [@Gralla] with some minor changes made to suit the initial conditions at infinity.

We expect this paclet to be used primarily for studying the observational characteristics of accretion disks of black holes, including dynamical processes in the disk, since the time delay of each geodesic is also computable in this paclet. The paclet can also be used for educational purposes, such as plotting null geodesics of the Kerr metric or illustrating the lensing effects, or for any other purpose where the computation of photon trajectories is necessary.

The paclet may be complemented with the `KerrGeodesics` Mathematica package [@BHPToolkit] which can be used to compute timelike geodesics. Among published codes that deal in more detail with the physics of accretion disks, we might mention [Tlusty]{.sc}, a Fortran77 code that can be used to model stellar atmospheres and accretion disks [@Hubeny].