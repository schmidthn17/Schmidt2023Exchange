# Schmidt2023Exchange
Primary code files that generate the data for the paper "Exchange between compartments regulates steady states and stochastic switching of a multisite phosphorylation network."

Generating data: TwoCompartmentStochasticSimulation
  Total particle(s) and their distribution can be altered in lines 22 - 39 of the above file
  Number of trajectors (nTraj), max time for the simulation (tMax), and the size of the time step where you collect data (deltaTime) are tunable
  Compartmental volume is tunable using parameters "vA" and "vB"
  Compartmental exchange rates for each species are tunable using lines 70 - 87 of the above file

Identifying stochastic switching events and the spread of switching events for each data set: 
