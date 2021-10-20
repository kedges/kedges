# Optimization on the Smallest Eigenvalue of the Grounded Laplacian Matrix via Edge Addition

Julia code for the article "Optimization on the Smallest Eigenvalue of the Grounded Laplacian Matrix via Edge Addition"


[edgecore.jl](./edgecore.jl) and [alg1.jl](./alg1.jl) contain Greedy and Fast and other baselines mentioned in the article, which are function exactgreedy() and function greedy() and so on.

Algorithm Vector mentioned in the article to calculate the smallest eigenvector can be found in [edgecore.jl](./edgecore.jl) and is the function eigcal()


[edgemain.jl](./edgemain.jl) is the main body of the algorithm and calls function exactgreedy() and function greedy() and other baselines to compare their effectiveness and efficiency


[soc-dolphins.mtx](./data/soc-dolphins.mtx) is an example of the data
