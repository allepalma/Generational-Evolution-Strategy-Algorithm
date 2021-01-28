from IOHexperimenter import IOH_function, IOH_logger, IOHexperimenter
import numpy as np
import sys

budget = 10000


'''Recombination functions'''

def global_discrete(a,lamb,n):
    #Initialize the vector that will host the offspring
    offspring = []
    #For as many times as many offspring we have
    for _ in range(lamb):
        off = []
        for i in range(n):
            #Recombine the individuals
            ind = np.random.choice(range(mu))
            off.append(a[ind][0][i])
        #Give to the offspring the sigma from the last randomly picked individual.
        sigma_rec = a[ind][1]
        #Add the recombined lambda and recombined individual 
        offspring.append((np.array(off),sigma_rec))
    return offspring

def global_inter(a, lamb):
    #Initialize the average individual with zeros
    off = np.zeros(d)
    #Initialize the average sigma with zeros
    sigma_rec = 0
    for ind in a:
        off += ind[0]
        sigma_rec += ind[1]
    couple = (off/mu, sigma_rec/mu)
    return [couple for _ in range(lamb)]

def pairs_discrete(a,lamb, n):
    #The total offspring of size mu
    offspring = []
    #Set up the indexes
    indexes = np.arange(mu)
    for _ in range(lamb):
        off = []
        #Pick randomly indexes of two parents
        parent1 = np.random.choice(indexes)
        parent2 = np.random.choice(np.delete(indexes,parent1))
        pair = (parent1,parent2)
        for i in range(n):
            ind = np.random.choice([0,1])
            off.append(a[pair[ind]][0][i])
        #Again, give to it the sigma of the last random parent picked 
        sigma_rec = a[pair[ind]][1] 
        offspring.append((np.array(off),sigma_rec))
    return offspring


def pairs_inter(a, lamb):
    #Initialize offspring
    off = []
    indexes = np.arange(len(a))
    for _ in range(lamb):
        parent1 = np.random.choice(indexes)
        parent2 = np.random.choice(np.delete(indexes,parent1))
        off.append(((a[parent1][0]+a[parent2][0])/2, (a[parent1][1]+a[parent2][1])/2))
    return off


'''Mutation'''

#Functions to mutate sigma
def mutate_one_sigma(sigma,tao0):
    return sigma*np.exp(np.random.normal(0,tao0, len(sigma))

#The individual sigma strategy
def mutation_one(a, tao):
    off = []
    for i in range(len(a)):
        sigma_prime = mutate_one_sigma(a[i][1], tao)
        ind_prime = a[i][0] + np.random.normal(0,sigma_prime)
        off.append((ind_prime, sigma_prime))
    return off
        

def optimiz(problem):
    #Set the total number of variables (dimensionality)
    n = problem.number_of_variables

    #At the beginning the optimum is the maximum system integer
    fopt = -sys.maxsize-1

    #Initialize the population and its standard deviations
    x = [np.random.rand(n) * 10 - 5 for _ in range(mu)]
    sigma = np.repeat(0.5, repeats = mu)
    a = list(zip(x,sigma))

    
    #Initialize the tao variable used for the self-adaptation 
    tao0 = 1/(np.sqrt(n))


    #Start the central loop
    while not problem.final_target_hit and problem.evaluations < budget:
        
        #Recombine
        if mu > 1:
            if rec == 'glob':
                if rec_type == 'discrete':
                    a_prime = global_discrete(a, lamb, n)
                else:
                    a_prime = global_inter(a, lamb)
            elif rec == 'pairs':
                if rec_type == 'discrete':
                    a_prime = pairs_discrete(a, lamb, n)
                else:
                    a_prime = pairs_inter(a, lamb)
        if mu == 1:
            a_prime = [a[0] for _ in range(lamb)]
            

        #Mutate
        a_prime_prime = mutation_one(a_prime, tao0)

        
        #Select 
        if sel == 'mu+lamb':
            a_tot = a+a_prime_prime

        else:
            a_tot = a_prime_prime
    
            
        f = {i:problem(a_tot[i][0]) for i in range(len(a_tot))}
        f_sorted = sorted(f.items(), key=lambda item: item[1])
        a = [a_tot[i[0]] for i in f_sorted[:mu]]
        fopt = f_sorted[0][1]
    return a, fopt





if __name__ == '__main__':

    ## Declarations of Ids, instances, and dimensions that the problems to be tested.
    problem_id = range(1,25)
    instance_id = range(1,26)
    dimension = [2, 5, 20]

    mu = 2
    lamb = 3

    #rec between glob and pairs
    rec = 'glob'
    #rec_type between interm and global
    rec_type = 'inter'
    #sel mu+lamb or mu/lambd
    sel = 'mu,lamb'
    
    
    ## Declariation of IOHprofiler_csv_logger.
    logger = IOH_logger("./", "one_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,"one_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,
                        "one_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type)

    for p_id in problem_id :
        for d in dimension :
            for i_id in instance_id:
                ## Getting the problem with corresponding id,dimension, and instance.
                f = IOH_function(p_id, d, i_id, suite="BBOB")
                f.add_logger(logger)
                xopt, fopt = optimiz(f)
    logger.clear_logger()
