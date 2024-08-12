# REEFMOD.7.0_GBR

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR).

Version 7.0 enables future projections (2008-2100) under a suite of CMIP-6 climate models for five greenhouse gas emission scenarios (SSP1-1.9, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5) with the simulation of mechanisms of coral adaptation (not optional). This version has been used for producing the 2024 counterfactual scenarios of the Reef Restoration and Adaptation Program (RRAP: https://gbrrestoration.org/).

Either CMIP-5 or CMIP-6 models can be selected from the front script MAIN_REEFMOD_GBR.m (but see REEFMOD.6.8_GBR/MAIN_REEFMOD_GBR.m for using CMIP-5 models as input). The script is designed to run the model under one climate change scenario (ie, one CMIP-5 or CMIP-6 model with the available scenario of carbon emission RCP/SSP) chosen by the user. The number of repeat simulations, ie, the stochastic simulation of the same warming scenario under different projections of cyclones and other random events (including the initialisation of coral cover and CoTS density, the magnitude of coral mortality, the forcing scheme of water quality,...) can be set with 'NB_SIMULATIONS' (set to 20 in RRAP). Because the runtime of one model run is about 2 hours, use of HPC resources is recommended.

The model reconstructs coral trajectories across the GBR between 2008-2024 and forecasts possible coral futures (2024-2100) based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (mass coral bleaching) and the simulated population dynamics of the coral-eating crown-of-thorns starfish (CoTS). Options for simulating management interventions on any given reef include:
- the outplanting of corals of specified species group, size and heat tolerance at a specified density
- the enrichment of coral larvae of specified species group with a specified density
- the reduction of heat stress through Solar Radiation Management (fogging)
- the consolidation of lose coral rubble to increase survival of coral recruits

This version also simulates the CoTS control program in space and time from 2019 onwards, with a specific number of boats (default: 5 boats) and the list of priority reefs as currently (2021) defined by GBRMPA.

Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)

Citation: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2022. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494
https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494
