---
title: 'DYCOVE: A Python package for coupling dynamic vegetion processes with hydro-morphodynamic models'
tags:
  - Python
  - vegetation
  - eco-morphodynamics
  - salt marshes
  - numerical modeling
authors:
  - name: Nelson Tull
    orcid: 0000-0002-6690-1106
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Muriel Brückner
    orcid: 0000-0002-7954-9586
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Civil and Environmental Engineering, Louisiana State University
   index: 1
 - name: Center for Computation & Technology, Louisiana State University
   index: 2
date: 14 December 2025
bibliography: paper.bib
---

# Summary

Vegetation growth in coastal environments plays an important role in shaping coastal morphology `[@jones1994; @kirwan2016; @kleinhans2018; @mariotti2010; @schwarz2018; @temmerman2005; @temmerman2007]`. Hydrodynamic and morphodynamic (numerical) models are used widely for understanding the processes that impact coastal systems, and they inform management strategies for coastal protection and ecosystem preservation. These models typically have options to incorporate constant vegetation effects via a roughness parameter, like Manning's $n$ or Chézy-value, where taller and denser vegetation is associated with higher roughness that alters flow velocities and sediment transport across landscapes. However, vegetation in these environments is not static, and their attributes are not constant: the colonization, growth, and mortality of a vegetation species is a dynamic process that depends on environmental conditions and species interactions. As a result, dynamic vegetation representations in numerical models are necessary to determine their feedbacks with local hydrodynamics and morphodynamics. This work presents an open-source Python modeling framework DYCOVE (DYnamic COastal VEgetation) that is dynamically coupled with the commonly used numerical models Delft3D Flexible Mesh (DFM) and ANUGA, with capabilities to be expanded to other models. DYCOVE allows us to represent life-cycle dynamics of multiple vegetation species in coastal environments through detailed vegetation processes that are functions of environmental conditions, resulting in spatial and temporal updates of friction effects in the numerical model. The model reproduces realistic vegetation patterns and roughness distributions and allows us to test how varying degrees of vegetation complexity affect hydro-morphodynamics in coastal systems. This project is designed to prioritize straightforward modeling of coastal vegetation dynamics, accessibility for all user levels, and smooth implementation of additional ecological and biological processes relevant to the study of biophysical interactions.


# Statement of Need

DYCOVE's vegetation logic originated as a MATLAB model `[@bruckner2019]` coupled with Delft3D-FLOW `[@lesser2004]`. In that study, the dynamic vegetation approach was shown to reproduce observed vegetation distributions of estuarine salt marshes. This original model was  as one of the first models that incorporated detailed ecological processes and biophysical feedbacks in large-scale hydro-morphodynamic numerical models, but it relieson proprietary software (MATLAB) and is limited in its coupling capability with other numerical models. DYCOVE is an advanced Python-implementation of the original model, which is now object-oriented, completely open-source, and has been reorganized and documented in a way that will allow easy implementation ofdynamic vegetation processes into coastal models. Currently, DYCOVEcan couple with the widely-used and modernized Delft3D Flexible Mesh Suite `[@delft3d]` while also offering a fully open-source avenue with the ANUGA model `[@anuga]`. Furthermore, DYCOVE is written in such a way that will allow for straightforward integration of other hydro-morphodynamic models and modules in the future, for example, pyDeltaRCM `[@moodie2021]`.


# Background

DYCOVE allows users to incorporate coastal vegetation in their numerical models through a dynamic coupling framework. I would add some detail here about the coupling: feedbacks at dense time-scales etc. Then next the equations that represent the dynamic vegetation (growth, mortality etc). Then the part on how to do it:
set of vegetation species input parameters must be provided in the form of a `.json` file that describes species attributes such as initial and maximum stem height and diameter, number of life stages, colonization window, and more. Simulations are run by calling a model engine class (depending on which numerical model is being used) that implements the time stepping methods specific to that model. At regular intervals, typically one tidal cycle, hydrodynamic and morphodynamic statistics from the preceding interval are used to compute changes in vegetation state (colonization and mortality). Vegetation stem height, stem diameter, and root length growth are determined purely based on species characteristics defined in the input `.json` file. These parameters, along with the species' stem density, are passed back to the numerical model via the Baptist formulation `[@baptist2007]`, which calculates a Chezy roughness value:

$$C_v=\frac{1}{\sqrt{(C_b^{-2} + \frac{C_d m D}{2 g} \min(h, h_v))}} + 
\frac{\sqrt{g}}{\kappa} \ln{\frac{\max(h, h_v)}{h_v}}$$

where $C_v$ is the vegetated Chezy coefficient, $C_b$ is the bare-earth Chezy coefficient,  $C_d$ is the vegetation drag coefficient, $m$ is stem density [$\text{m}^{-2}$], $D$ is stem diameter [m], $h$ is flow depth [m], $h_v$ and stem height [m], $g$ is gravitational acceleration [$\text{m}/\text{s}^2$], and $\kappa$ is the von Kármán constant.

DYCOVE simulations run on the idea that vegetation and morphological changes both occur on longer time scales than hydrodynamic changes. Therefore, we "accelerate" vegetation processes in a manner similar to the Morphological Acceleration Factor (MORFAC) of DFM `[@ranasinghe2011]`. In fact, for DYCOVE-DFM simulations with morphology active, we set the Ecological Acceleration Factor (``ecofac``) equal to the ``morfac`` of the DFM model.

To reconcile the hydrodynamic and ecological time scales, we run simulations by defining the following parameters:

- ``ecofac``: Ecological Acceleration Factor, a multiplier representing the ratio of ecological time to hydrodynamic time.
- ``veg_interval``: Vegetation coupling interval, in hydrodynamic seconds, defining how often changes in vegetation are computed and passed back to the hydrodynamic model.
- ``n_ets``: Number of Ecological Time Steps, or vegetation coupling intervals, per ecological year.

Because DYCOVE is most useful for modeling vegetation in the intertidal zone, it is convenient to structure simulations around full tidal cycles. For your typical semi-diurnal tide, one cycle is about 12 hours (43200 seconds), rounding down for simplicity. These variables relate to each other in the following way:

$$\text{ecofac} \approx \frac{365 \times 86400}{\text{veg interval} \times \text{n}_{\text{ets}}}$$


# Easy Implementation


# Acknowledgements

We acknowledge the Center for Computation and Technology for funding this work.
We also thank Vindhyawasini Prasad and Anderson Amaya Saldarriaga for their help in testing the package and providing ideas.


# References