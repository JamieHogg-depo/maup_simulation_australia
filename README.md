# maup_simulation_australia

This repository contains code for the analyses described in the manuscript "The effect of the modifiable areal unit problem on ecological model inference: A simulation study for disease mapping in Australia" by James Hogg, Aiden Price, Conor Hassan, Shovanur Haque, Farzana Jahan, Wala Areed, Jessica Cameron, Helen Thompson, Susanna Cramb

## Notes

The simulation experiment was carried out using a high performance computer, by running `run_sim.sh`. This file creates many `.rds` files which are then combined using `combine_results.sub`. The code to create some plots can be found in `getResultsPlots.R`

The remainder of the files are available [here](https://drive.google.com/drive/folders/12vM5dS96Bu8h_juTgOVrt9LP87GRRG5N?usp=sharing)
- phi
- sa1_shape
- zones
- rds

The three folders should be added into `r_src` prior to running the code. 
