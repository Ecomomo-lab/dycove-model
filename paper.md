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

Vegetation growth in coastal environments plays an important role in shaping coastal morphology [@jones1994; @kirwan2016; @kleinhans2018; @mariotti2010; @schwarz2018; @temmerman2005; @temmerman2007]. Hydrodynamic and morphodynamic (numerical) models are used widely for understanding the processes that impact coastal systems, and they inform management strategies for coastal protection and ecosystem preservation. In these models, vegetation effects are incorporated via friction, with increasing roughness effects for taller and denser vegetation that alters flow velocities and sediment transport. However, roughness is commonly assumed as being constant, a key limitation as vegetation is not static in space and time: processes such as colonization, growth, mortality and species interactions determine dynamic patterns of biomass that depend on environmental conditions. As a result, dynamic vegetation representations are necessary to determine their interactions with local hydrodynamics and morphodynamics in numerical models. This work presents an open-source Python modeling framework DYCOVE (DYnamic COastal VEgetation) that is dynamically coupled with the numerical models Delft3D Flexible Mesh (DFM) and ANUGA, with capabilities to be expanded to other models. DYCOVE represents life-cycle dynamics of multiple vegetation species in coastal environments through detailed ecological processes that are functions of environmental conditions, resulting in spatial and temporal updates of friction effects in the numerical model. The model reproduces realistic vegetation patterns and roughness distributions, a prerequisite to test howf vegetation complexity affects hydro-morphodynamics in coastal systems. DYCOVE is designed to prioritize straightforward modeling of coastal vegetation dynamics, accessibility for all user levels, and smooth implementation of additional ecological and biological processes relevant to the study of biophysical interactions.
 
# Statement of Need

As the scientific community has become increasingly aware of the relevance of biophysical feedbacks for coastal change, there has been a growing focus on understanding the role of vegetation in coastal protection and ecosystem function (e.g., salt marshes) and hence their incorporation in predictive models [@bruckner2019; @gourgue2022; @jeanlouis2025, @wright2018]. With this interest in mind, we present a dynamic vegetation model that is coupled with commonly used hydro-morphodynamic modeling packages.

DYCOVE's vegetation logic originated from a MATLAB model [@bruckner2019] coupled with Delft3D-FLOW [@lesser2004] that successfully reproduced observed vegetation distributions of estuarine salt marshes without the need of calibration. This original model was one of the first to incorporate detailed ecological processes and biophysical feedbacks in large-scale hydro-morphodynamic numerical models, but it relies on proprietary software (MATLAB) and is limited in its coupling capability with other numerical models. DYCOVE is an advanced Python-implementation of the original model, which is now object-oriented, completely open-source, and has been reorganized and documented in a way that will allow easy implementation of dynamic vegetation processes across coastal models. 
 
# Background
DYCOVE incorporates coastal vegetation in physics-based numerical models through a dynamic coupling framework, where vegetation is distributed based on habitat parameters and provides spatially-varying roughness. The response of the physics-based model to this roughness leads to alteration of the environmental conditions, i.e., inundation, flow velocities, etc., that define habitat and vegetation is updated during the next coupling, resulting in a biophysical feedback-loop.

Vegetation dynamics are simulated through rules for colonization, growth, senescence, and mortality. Colonization occurs during a defined time step when wetted cells receive a user-defined vegetation fraction, after which vegetation ages and grows linearly up to specified maximum dimensions, with optional winter stem height to represent aboveground biomass mortality. Senescence is modeled via a maximum age, supporting both annual and perennial vegetation, and optional life stages allow vegetation properties to change as plants mature. Mortality is a dose-effect relationship where vegetation fractions are gradually removed if habitat conditions are not optimal. The removed vegetation fraction depends linearly on the strength of the mortality pressure.
 
The vegetation module is coupled to the numerical model at an interval that we define as the Ecological Time Step (ETS), which is typically set to one tidal cycle (ref bruckner et al). At the end of each ETS, the numerical model provides time-varying quantity information (inundation fraction, flow velocities, and bed level changes) that drives vegetation colonization and mortality. Vegetation attributes are passed back to the numerical model via the Baptist formulation [@baptist2007], which calculates a Chézy roughness value:
 
$$C_v=\frac{1}{\sqrt{(C_b^{-2} + \frac{C_d m D}{2 g} \min(h, h_v))}} +
\frac{\sqrt{g}}{\kappa} \ln{\frac{\max(h, h_v)}{h_v}}$$
 
where $C_v$ is the vegetated Chézy coefficient [$\text{m}^(½) \text{s}^(-1)$], $C_b$ is the Chézy coefficient for alluvial bed roughness [$\text{m}^(½) \text{s}^(-1)$], $C_d$ is the vegetation drag coefficient, $m$ is stem density [$\text{m}^{-2}$], $D$ is stem diameter [m], $h$ is flow depth [m], $h_v$ and stem height [m], $g$ is gravitational acceleration [$\text{m}/\text{s}^2$], and $\kappa$ is the von Kármán constant. To achieve meaningful changes in vegetation over the course of a simulation, we "accelerate" vegetation processes by a constant scaling factor, using the same logic as the MORFAC parameter in Delft3D [@ranasinghe2011].
 
# Model Implementation

DYCOVE currently supports coupling to DFM and ANUGA, and therefore relies on installation of one (or both) of those models along with the DYCOVE source code. As such, DYCOVE is geared toward users who are familiar with one of these models or have an existing model to which they would like to add dynamic vegetation processes.

Users who have little to no experience with vegetation modeling, DFM, or ANUGA can take advantage of the DYCOVE-ANUGA workflow that can be installed directly from GitHub and have one of the provided working examples running within a matter of minutes. On the other hand, support for DFM carries the benefit of familiarity to many users in the coastal modeling community, while also having a morphodynamic component that is an important aspect of ecological modeling in coastal environments [@bruckner2019].

Regardless of the underlying numerical model, DYCOVE is implemented with just a few lines of code, without any need to modify the model setup users may already have with DFM or ANUGA. The modularity of DYCOVE’s underlying numerical model classes will allow for straightforward integration of other hydro-morphodynamic models in the future. For example, pyDeltaRCM [@moodie2021] is an open-source and well-documented option with a morphology component that could be added to DYCOVE without a major code refactoring.


# Acknowledgements

We acknowledge the Center for Computation and Technology for funding this work.
We also thank Vindhyawasini Prasad and Anderson Amaya Saldarriaga for their help in testing the package and providing ideas.


# References