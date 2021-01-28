# Generational-Evolution-Strategy-Algorithm

This repository is dedicated to the implementation of a canonical generational Evolutionary Strategy algorithm for the optimization of complex black-box objective functions defined over a multi-dimensional real space. Candidate solutions are generated as n-dimensional vectors and treated through standard genetic operations of recombination, mutation and selection. The benefits and drawbacks of different approaches to set up such operators can be evaluated in different contexts. We applied our algorithm to 24 noise-free real-parameter single-objective benchmark functions (BOBB suite) with different characteristics and properties implemented in the python IOHexperimenter package for heuristic optimization. To further explore the potentiality of our algorithm, we examined the performance of the objective function ensemble on three different space dimensions: 2-dimensional, 5-dimensional and 20-dimensional. Even if carefully tailored to our setting, our approach strives towards re-usability and our functions can be easily adapted to multiple scripts and contexts.

## The implementation 
Three different versions of the Evolution Strategy algorithm are being proposed exploiting three distinct candidate solution representation and genetic operators: the one-sigma strategy, the individual sigma strategy and the individual sigma strategy with correlated dimensions. The algorithms are run across three dimensionality and aiming at 25 different target values. The scope of the optimization algorithm is to find a candidate solution vector such when the objetive functions of the suite are evaluated on it they output which is as close a as possible to a pre-defined target.

### The representation
The individuals in the population are randomly generated solutions (with all dimensions embedded in the interval [-5,5]). Each individual is made of a vector of real numbers associated with a step-size vector indicated with the Greek letter sigma. The step size acts as standard deviation of a normal distribution centred at 0 from which random variates are drawn at each step to update the multi-dimensional vector of solutions. In the one-sigma representation, a single step-size value is used to control all dimensions of the vector of solutions. In the individual sigma strategy, the step-size is a vector with the same dimensionality as the vector of solutions, such that the ith element of the step-size vector parametrizes the normal random variable used to update the ith dimension of the condidate solution. The last approach, in addition to a step-size vector, binds each solution with a vector of <img src="https://render.githubusercontent.com/render/math?math=\frac{n \times (n-1)}{2} "> rotational angles between each pair of dimension. The mutation magnitude and direction in this last scenario is obtained as a vector of variates from a multi-dimensional normal distribution parametrized by a covariance matrix.


### The recombinations
The two types of recombination we implemented are intermediary recombination and discrete recombination. Given a pool of parent solutions to recombine for the production of new individuals in light of the mutation step, the intermediary recombination generates offsprings by averaging the vectors of parent individuals. Conversely, the discrete recombination produces new solutions by inserting at a certain position the value presented by a random parent of the parental pool at that specific position. Now it remains to define how to select the parent pool. In this case, two alternative strategies were again enacted for such a process to be carried out. The first (and faster) approach was
the one of global recombination, which implements the rules of intermediary and discrete recombination to the parent solutions made by all the individuals in the parent population. Alternatively, we also experimented with canonical recombination between pairs of parents. In this latter case, instead of considering all individuals we draw couples of parents and use them for recombination.

### The mutation
The mutation operator regards both the solution vector and the model parameters (step-size/rotational angles). For each iteration of the optimization process, first the sigma vector is mutated through a random variate drawn from a log-normal distribution (depending on a single learning rate parameter in the one-sigma scenario and both a global and local learning rate in the two more complex settings). Subsequently, the candidate solution is updated by adding it to a random vector drawn from a random normal distribution parametrized by the previously mutated step-size vector/value and centred at 0. In case we implement correlated mutations too, before modifying a candidate solution during a run, we first mutate the rotational angle vector too through a random vector from a normal distribution parametrized by a constant and centred at 0. In this last situation, before performing any change to the central vector of solutions we apply pairwise rotations between each pair of dimensions of the random mutation vector via the product of <img src="https://render.githubusercontent.com/render/math?math=\frac{n \times (n-1)}{2} "> rotational matrices.

### The selection operator
 







