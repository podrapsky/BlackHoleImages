---
title: 'BlackHoleImages: A Wolfram Mathematica paclet for analytical black hole imaging'
tags:
  - Accretion Disks
  - Analytical Methods
  - Black Holes
  - Relativistic Astrophysics
  - Wolfram Mathematica
authors:
  - name: David Podrápský
    orcid: 0009-0000-7558-2509
    affiliation: 1
  - name: Vojtěch Witzany
    orcid: https://orcid.org/0000-0002-9209-5355
    affiliation: 1
affiliations:
  - name: Institute of Theoretical Physics, Faculty of Mathematics and Physics, Charles University, V Holešovičkách 2, 180 00 Praha 8, Czech Republic 
    index: 1
date: 18 February 2025
bibliography: paper.bib

---

# Summary

Black holes are astrophysical objects with gravitational fields so extreme that they significantly affect even particles moving at the speed of light in their vicinities. Consequently, the presence of a black hole can lead to strong lensing effects on the images emerging from their surroundings. In the geometrical optics limit, the paths of particles and fields moving at the speed of light follow null geodesics in the spacetime geometry. Additionally, the passage of light through the gravitational field leads to changes in the frequency and intensity of the light, which is particularly important to understand the observational imprints of accretion disks forming around black holes. 

`BlackHoleImages` is a Wolfram Mathematica paclet designed to compute the abovementioned effects in the field of isolated rotating black holes in Einstein's General Relativity, which are described by the so-called Kerr family of spacetimes. The paclet can be used for computing: 
- photon orbits (null geodesics) in the field of black holes,
- generating images of stellar backgrounds distorted by black holes, and 
- generating images of accretion disks near black holes.

The paclet contains three packages, `KerrNullGeodesics`, `AlphaDiskModel`, and `KerrImages`. `KerrNullGeodesics` implements the analytical solutions for null geodesics in Kerr space-times as provided by Gralla and Lupsasca [@Gralla].  The `AlphaDiskModel` then implements the standard $\alpha$-disk model for the accretion disk near a black hole as proposed by Shakura and Sunyaev [@Shakura], specifically its relativistic incarnation by Novikov & Thorne [@NovikovThorne] and Page & Thorne [@PageThorne]. Finally, the `KerrImages` package uses the previous two packages to provide functions that compute the disk images and the lensed images of faraway backgrounds.

We expect this paclet to be used primarily for investigating the observational characteristics of accretion disks of black holes, including dynamical processes in the disk, since the time delay of each geodesic is also computable in this paclet. The paclet can also be used for educational purposes, such as plotting null geodesics in the Kerr metric or illustrating the lensing effects, or for any other purpose where the computation of photon trajectories is necessary.

# Statement of need

Black holes are one of the most intriguing objects in modern physics. They were on of the earliest predictions of Einstein's general relativity and started to be indirectly through observations in various electromagnetic bands starting from the 1960s [@EMBHRev]. Since most observations were and still are unable to resolve the spatial scales of the so-called event horizon of the black hole, we need to understand the radiation that emerges from the interactions of the black hole with its surroundings and the emergent spectrum. This was provided by the seminal $\alpha$ accretion disk model of Shakura and Sunyaev [@Shakura] and Novikov & Thorne [@NovikovThorne]. As a result, the community of astronomers and astrophysicists identified the first candidates for black holes in the Universe.

Currently, the evidence for dark compact objects consistent with black holes is overwhelming. Throughout the 1990s, two teams of astronomers followed stars orbiting very close to our Galactic center and deduced on the presence of a dark compact object with the mass of $4\cdot 10^6 M_\odot$ ($M_\odot = $ solar mass) [@Genzel,@Ghez]; the leaders of the teams later shared the 2020 Nobel prize in physics [@NobelGenzelGhez]. In 2015, the LIGO collaboration detected the first merger of black holes by gravitational waves [@GWBH] and more than 180 events were detected to date [@GWOSC]. In 2017, the Event Horizon Telescope colaboration reconstructed the image of the accretion flow in the immediate vicinity of the massive black hole in the center of the M87 galaxy by the methods of very long baseline interferometry [@EHTM87]. They followed up this result in 2020 by presenting a similar image of the accretion flow from the center of our very own Galaxy the Milky Way  [@EHTSgrA].

The images of the massive black holes from M87 and the Milky Way present a dark feature in their centers known as the ``shadow'' of the black hole. It corresponds to a region of the sky where the paths of the photons are bent so much that they end behind the black hole horizon, the point from beyond which not even light can escape. While this qualitative picture is correct, the precise quantitative science inferred from the observations of the Event Horizon Telescope requires the detailed modeling of the accretion flow, its radiation, and the emerging photon orbits [@EHTmodelling]. This is one of the important use cases our paclet is aiming for. 

Another important application of the paclet is the interpretation of timing and spectral features of unresolved accreting black holes as observed in the burgeoning field of time-domain astronomy [@TDAMM]. Transients such as quasi-periodic eruptions, oscillations, or outflows emerging from systems containing accreting black holes are detected at an increasing rate. They often contain subtle timing features that require precise modelling of the accretion flows and other physical components, but possibly also the lensing effects by the black holes. Our paclet can be easily extended to compute the lightcurves and/or images emerging from these systems.

Last but not least, the topic of black holes captures the imagination of the public. In the 2014 movie *Interstellar* directed by Christopher Nolan, the general movie audience was able to observe a faithful simulation of a black hole image for the first time [@InterstellarModeling]. Our paclet provides an intermediate step for relativity students and other interested parties to be able to disect how such an image is made through a user-friendly Wolfram Mathematica paclet.

On the more technical side, this is also one of the first open source releases of an implementation of the *analytical* solution for null geodesics in Kerr space-time (the only other implemenation we know of is [@NullGeodesicsBranch]). Other implementations typically use numerical differential equation solvers (e.g. the EinsteinPy implementation [@EinsteinPyKerrGeo]). The arbitrary precision capabilities of *Wolfram Mathematica* thus allow to use this paclet as an exact validation and reference solution for such machine precision codes. Using the symbolic calculus capabilities of Mathematica provided a good platform to implement the relatively complicated analytical solution as a starting point. The paclet could now also be used as a reference when porting the solutions to more optimized lower-level implementations. For example, the Mathematica code `KerrGeodesics` [] for time-like geodesics (trajectories of massive test particles) in Kerr space-time was recently rewritten in Python as `KerrGeoPy` [], our paclet could follow a similar route.  
