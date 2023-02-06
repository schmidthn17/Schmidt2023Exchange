# Schmidt2023Exchange
Primary code files that generate the data for the paper "Exchange between compartments regulates steady states and stochastic switching of a multisite phosphorylation network."

Generating data: TwoCompartmentStochasticSimulation.m
  Total particle(s) and their distribution can be altered in lines 22 - 39 of the above file
  Number of trajectors (nTraj), max time for the simulation (tMax), and the size of the time step where you collect data (deltaTime) are tunable
  Compartmental volume is tunable using parameters "vA" and "vB"
  Compartmental exchange rates for each species are tunable using lines 70 - 87 of the above file
  Do NOT change parameters k1 - k6 as these are experimentally measured values 

Identifying stochastic switching events: IdentifyingStochasticSwitches.m 
  This file generates an mp4 of all trajectories with switches identified on each plot for both compartments. Calculations are done to provide the min, mean, and max number of switching events.
  Line 36 holds the file name for your mp4 file and should be edited for each data set to avoid overwriting your movie file. 
  No other edits to this file are necessary to run if you're following the parameterization from the paper.
