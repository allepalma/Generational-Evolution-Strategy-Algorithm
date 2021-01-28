from IOHexperimenter import IOH_function, IOH_logger, IOHexperimenter
import numpy as np
import sys
import time
import itertools
import copy

budget = 10000


'''Recombination functions'''

def pairs_discrete(x, sigma,lamb):
    #Initialize the recombinant matrices as filled with zeros
    sigma_rec = np.zeros((lamb,d))
    x_rec = np.zeros((lamb,d))
    #Fix the indexes of the rows of the matrix
    indexes = np.arange(len(x))
    for i in range(lamb):
        parent1 = (np.random.choice(indexes))
        parent2 =  np.random.choice(np.delete(indexes,parent1))
        parents = (parent1, parent2)
        for j in range(d):
            #Fetch one parent to fill up the recombinant rows.
            parent = parents[np.random.choice([0,1])]
            x_rec[i,j] = x[parent, j]
            sigma_rec[i,j] = sigma[parent, j]
    return x_rec, sigma_rec
                
                        
def global_discrete(x, sigma,lamb):
    #Initialize the recombinant matrix as filled with zeros
    sigma_rec = np.zeros((lamb, d))
    x_rec = np.zeros((lamb, d))
    #Cycle across the columns
    for i in range(x.shape[1]):
        sigma_rec[:,i] = np.random.choice(sigma[:,i], lamb)
        x_rec[:,i] = np.random.choice(x[:,i], lamb)
    return x_rec, sigma_rec

def global_inter(x, sigma,lamb):
    #Initialize the average individual with zeros
    x_rec = np.tile(np.mean(x, axis=0),(lamb,1))
    sigma_rec = np.tile(np.mean(sigma, axis=0), (lamb,1))
    return x_rec, sigma_rec


def pairs_inter(x, sigma, lamb):
    #Initialize offspring
    x_rec = np.zeros((lamb, d))
    sigma_rec = np.zeros((lamb, d))
    indexes = np.arange(len(x))
    #Set up the indexes of the rows of the matrix
    for i in range(lamb):
        parent1 = np.random.choice(indexes)
        parent2 = np.random.choice(np.delete(indexes,parent1))
        x_rec[i,:] = (x[parent1,:]+x[parent2,:])/2
        sigma_rec[i,:] = (sigma[parent1,:]+sigma[parent2,:])/2           
    return x_rec, sigma_rec


'''Mutation'''

#Functions to mutate sigma
def mutate_ind_sigma(sigma,tao_prime, tao):
    sigma_prime = sigma*np.exp(np.random.normal(0,tao_prime,lamb).reshape((-1,1))+np.random.normal(0,tao, sigma.shape))
    return sigma_prime


###The individual sigma strategy
def mutation_ind(x, sigma, tao_prime, tao):
    sigma_prime = mutate_ind_sigma(sigma, tao_prime, tao)
    x_prime = x+np.random.normal(0,sigma_prime)
    return x_prime, sigma_prime


def optimiz(problem):
    n = problem.number_of_variables
    
    fopt = -sys.maxsize-1

    x = np.array([np.random.rand(n) * 10 - 5 for _ in range(mu)])
    sigma = np.array([np.repeat(0.5, repeats = d) for _ in range(mu)])
  
    #Initialize the tao variable used for the self-adaptation 
    tao_prime = 1/(np.sqrt(2*n))
    tao = 1/(np.sqrt(2*np.sqrt(n)))

    #Start the central loop
    while not problem.final_target_hit and problem.evaluations < budget:
        #Recombine
        if mu > 1:
            if rec == 'glob':
                if rec_type == 'discrete':
                    x_prime, sigma_prime = global_discrete(x, sigma, lamb)
                else:
                    x_prime, sigma_prime = global_inter(x, sigma, lamb)
            elif rec == 'pairs':
                if rec_type == 'discrete':
                    x_prime, sigma_prime = pairs_discrete(x, sigma, lamb)
                else:
                    x_prime, sigma_prime = pairs_inter(x, sigma, lamb)

        elif mu == 1:
            x_prime, sigma_prime = np.tile(x, (lamb, 1)), np.tile(sigma, (lamb, 1))

        #Mutate
        x_prime_prime, sigma_prime_prime = mutation_ind(x_prime, sigma_prime,  tao_prime, tao)

        #x_prime_prime, sigma_prime_prime = x_prime, sigma_prime
        #Select 
        if sel == 'mu+lamb':
            x_tot = np.vstack((x, x_prime_prime))
            sigma_tot = np.vstack((sigma, sigma_prime_prime))

        else:
            x_tot = x_prime_prime
            sigma_tot = sigma_prime_prime

            
        f = {i:problem(x_tot[i,:]) for i in range(len(x_tot))}
        f_sorted = sorted(f.items(), key=lambda item: item[1])
        indexes = [i[0] for i in f_sorted[:mu]]
        fopt = f_sorted[0][1]
        x = x_tot[indexes,:]
        sigma = sigma_tot[indexes,:]
    return x, fopt





if __name__ == '__main__':

    ## Declarations of Ids, instances, and dimensions that the problems to be tested.
    problem_id = range(1,25)
    instance_id = range(1,26)
    dimension = [2, 5, 20]

    mu = 2
    lamb = 3

    #rec between glob and pairs
    rec = 'pairs'
    #rec_type between interm and global
    rec_type = 'inter'
    #sel mu+lamb or mu/lambd
    sel = 'mu,lamb'
    
    
    ## Declariation of IOHprofiler_csv_logger
    logger = IOH_logger("./", "ind_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,"ind_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,
                        "ind_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type)

    for p_id in problem_id :
        for d in dimension :
            for i_id in instance_id:
                ## Getting the problem with corresponding id,dimension, and instance.
                f = IOH_function(p_id, d, i_id, suite="BBOB")
                f.add_logger(logger)
                xopt, fopt = optimiz(f)
    logger.clear_logger()
