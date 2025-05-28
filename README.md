# REEFMOD.7.0_GBR

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR).

Version 7.0 enables future projections (2008-2100) under a suite of CMIP-6 climate models for five greenhouse gas emission scenarios (SSP1-1.9, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5) with the simulation of mechanisms of coral adaptation (not optional). This version has been used for producing the 2024 counterfactual scenarios of the Reef Restoration and Adaptation Program (RRAP: https://gbrrestoration.org/).

The model reconstructs coral trajectories across the GBR between 2008-2024 and forecasts possible coral futures (2024-2100) based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (mass coral bleaching) and the simulated population dynamics of the coral-eating crown-of-thorns starfish (CoTS). Options for simulating management interventions on any given reef include:
- the outplanting of corals of specified species group, size and heat tolerance at a specified density
- the enrichment of coral larvae of specified species group with a specified density
- the reduction of heat stress through Solar Radiation Management (fogging)
- the consolidation of lose coral rubble to increase survival of coral recruits

This version also simulates the CoTS control program in space and time from 2019 onwards, with a specific number of boats (default: 5 boats) and the list of priority reefs as currently (2021) defined by GBRMPA.

## Instructions

The code is written in MATLAB (2023b or earlier versions).
To execute the model:
1. Download all the necessary scripts and folders.
2. Add them to your current MATLAB path.
3. In the Command Window, type:
    > run('MAIN_REEFMOD_GBR.m')

This will start the simulation. The current settings run a projection of one climate change scenario (ie, one CMIP-5 or CMIP-6 climate model under a specific scenario of carbon emission RCP/SSP - as specified by the user) for the period 2008-2100, with CoTS control as the only management intervention (counterfactual simulation).

The number of repeat simulations can be set with 'NB_SIMULATIONS' (currently set to 20). Simulations are then executed sequentially, each identified by the iterator "simul" (eg, from 1 to 20), which sets set a specific seed for the MATLAB random number generator, ensuring reproducibility of the results. Each simulation is stochastic, incorporating several randomised components, including the timing of future heat stress within each decade, the selection of a specific scenario of future cyclones, the initialisation of coral cover and Crown-of-Thorns Starfish (CoTS) density, the magnitude of coral mortality events, the forcing scheme of water quality. Because the runtime of one complete simulation (ie, from year 2008 to year 2100) is about 2 hours, the use of HPC resources is recommended. Shorter simulations can be obtained by setting a lower number of 6-month time steps ('NB_TIME_STEPS').

## Citation
Bozec, Y.-M., A. A. Adam, B. Arellano-Nava, A. K. Cresswell, V. Haller-Bull, T. Iwanaga, L. Lachs, S. A. Matthews, J. K. McWhorter, K. R. N. Anthony, S. A. Condie, P. R. Halloran, J. C. Ortiz, C. Riginos, and P. J. Mumby. 2025. A rapidly closing window for coral persistence under global warming. bioRxiv. https://www.biorxiv.org/content/10.1101/2025.01.23.634487v1.full

## Earlier model versions for the GBR
Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2022. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494
https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494

Castro-Sanguino, C., Y.-M. Bozec, S. A. Condie, C. S. Fletcher, K. Hock, C. Roelfsema, D. A. Westcott, and P. J. Mumby. 2023. Control efforts of crown‐of‐thorns starfish outbreaks to limit future coral decline across the Great Barrier Reef. Ecosphere 14:e4580. https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.4580

Mason, R. A., Y.-M. Bozec, and P. J. Mumby. 2023. Demographic resilience may sustain significant coral populations in a 2° C-warmer world. Global Change Biology 29:4152–4160. https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16741

## Contact
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)

