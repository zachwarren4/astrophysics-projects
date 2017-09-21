Final ASTR 3800 Project. Data files too large for github repository as were some of the results. 
The goal was to take several million points of data (dark matter particle locations) and craft halos
based on a radius which extended from one point. If another point was within radius r, it was added to the halo.
The project seemed simple but initial attempts resulted in extremely long runtimes. Eventually, a cKDTree, which is a multi
dimensional tree, was used to map all points based on proximity to one another. This resulted in a much more efficient way of 
creating halos. 