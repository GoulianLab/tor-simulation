# tor-simulation
Python source code accompanying Carey et al. "Regulated stochasticity in a bacterial signaling network permits tolerance to a rapid environmental change"

On run, requests the mean number of bursts of TorT and TorS expression per cell generation and the number of times to run the simulation. Each simulation run generates a YFP fluorescence value (arbitrary scale) for every descendent of a single parent cell after 10 doublings, using a simplifying model of the behavior of a PtorCAD-yfp reporter. The program output is a text file with the every YFP value of every 10th-generation cell for every simulation run concatenated into a single column.

Theory and description of the implementation of the simulation are described in the associated publication.
