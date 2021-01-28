from IOHexperimenter import IOH_function, IOH_logger, IOHexperimenter
import numpy as np
import sys
import time
import itertools

budget = 10000


'''Recombination functions'''

def pairs_discrete(x, sigma, alpha, lamb, n):
    #Initialize the recombinant matrices as filled with zeros
    sigma_rec = np.zeros((lamb,n))
    x_rec = np.zeros((lamb,n))
    alpha_rec = np.zeros((lamb,(n*(n-1))//2))
    #Fix the indexes of the rows of the matrix
    indexes = np.arange(len(x))
    for i in range(lamb):
        #Pick two parents at random
        parent1 = (np.random.choice(indexes))
        parent2 =  np.random.choice(np.delete(indexes,parent1))
        parents = (parent1, parent2)
        for j in range(n):
            #Fetch one parent to fill up the recombinant rows.
            parent = parents[np.random.choice([0,1])]
            x_rec[i,j] = x[parent, j]
            sigma_rec[i,j] = sigma[parent, j]
            if j < alpha.shape[1]:
                alpha_rec[i,j] = alpha[parent, j]
        #Keep on filling the recombinant alpha
        if j < alpha.shape[1]:
            for k in range(j, alpha.shape[1]):
                parent = parents[np.random.choice([0,1])]
                alpha_rec[i,k] = alpha[parent, k]
    return x_rec, sigma_rec, alpha_rec
                
                        
def global_discrete(x, sigma,alpha,lamb,n):
    #Initialize the recombinant matrix as filled with zeros
    sigma_rec = np.zeros((lamb, n))
    x_rec = np.zeros((lamb, n))
    alpha_rec = np.zeros((lamb, (n*(n-1))//2))
    #Cycle across the columns
    for i in range(x.shape[1]):
        #Pick lambda parents at random and use them.
        random_parents = np.random.choice(np.arange(len(x)), lamb)
        sigma_rec[:,i] = sigma[random_parents, i]
        x_rec[:,i] = x[random_parents, i]
        if i < alpha.shape[1]:
            alpha_rec[:,i] = alpha[random_parents, i]
    if i < alpha.shape[1]:
        for j in range(i, alpha.shape[1]):
            random_parents = np.random.choice(np.arange(len(x)), lamb)
            alpha_rec[:,j] = alpha[random_parents, j]
    return x_rec, sigma_rec, alpha_rec


def global_inter(x, sigma,alpha, lamb):
    #Initialize the average individual with zeros
    x_rec = np.tile(np.mean(x, axis=0),(lamb,1))
    sigma_rec = np.tile(np.mean(sigma, axis=0), (lamb,1))
    alpha_rec = np.tile(np.mean(alpha, axis=0), (lamb,1))
    return x_rec, sigma_rec, alpha_rec


def pairs_inter(x, sigma, alpha, lamb):
    #Initialize offspring matrices
    x_rec = np.zeros((lamb, d))
    sigma_rec = np.zeros((lamb, d))
    alpha_rec = np.zeros((lamb, d*(d-1)//2))
    #Set up the indexes of the rows of the matrix
    indexes = np.arange(len(x))
    for i in range(lamb):
        parent1 = np.random.choice(indexes)
        parent2 = np.random.choice(np.delete(indexes,parent1))
        x_rec[i,:] = (x[parent1,:]+x[parent2,:])/2
        sigma_rec[i,:] = (sigma[parent1,:]+sigma[parent2,:])/2
        alpha_rec[i,:] = (alpha[parent1,:]+alpha[parent2,:])/2
    return x_rec, sigma_rec, alpha_rec

'''Correlate sigma'''

def correlate_sigma(sigma,alpha):
    n = d
    #n_q is the index of the rotational angle used at each iteration
    n_q = ((n-1)*n)//2
    sigma_corr = np.random.normal(0, sigma)
    for k in range(1,n):
        n_1 = n-k-1
        n_2 = n-1
        for i in range(k):
            print(n_2,n_1)
            d_1 = sigma_corr[:,n_1]
            d_2 = sigma_corr[:,n_2]
            sigma_corr[:,n_2] = d_1*np.sin(alpha[:,n_q-1]) + d_2*np.cos(alpha[:,n_q-1])
            sigma_corr[:,n_1] = d_1*np.cos(alpha[:,n_q-1]) - d_2*np.sin(alpha[:,n_q-1])
            n_2 = n_2-1
            n_q = n_q -1
    return sigma_corr
        

        
'''Mutation'''

#Functions to mutate sigma
def mutate_ind_sigma(sigma,tao_prime, tao):
    sigma_prime = sigma*np.exp(np.random.normal(0,tao_prime,lamb).reshape((-1,1))+np.random.normal(0, tao, sigma.shape))
    return sigma_prime


#The individual sigma strategy
def mutation_ind(x, sigma, alpha, beta, tao_prime, tao):
    #Implement the mutation of the individuals and parameters
    x_prime = x
    sigma_prime = mutate_ind_sigma(sigma, tao_prime, tao)
    alpha_prime = alpha+np.random.normal(0,beta,alpha.shape)
    sigma_prime_corr = correlate_sigma(sigma_prime, alpha_prime)
    x_prime = x_prime+sigma_prime_corr
    return x_prime, sigma_prime, alpha_prime

'''The algorithm implementing the reproduction cycle''' 

#Optimization algorithm
def optimiz(problem):
    #Initialize the number of variables
    n = problem.number_of_variables
    #Set the initial optimum
    fopt = -sys.maxsize-1
    #x is the individual, sigma the step size and alpha contains the rotation angles
    x = np.array([np.random.rand(n) * 10 - 5 for _ in range(mu)])
    sigma = np.array([np.repeat(0.5, repeats = d) for _ in range(mu)])
    alpha = np.array([np.repeat(0, repeats = (d*(d-1))//2) for _ in range(mu)])
    

      
    #Initialize the tao variable used for the self-adaptation 
    tao_prime = 1/(np.sqrt(2*n))
    tao = 1/(np.sqrt(2*np.sqrt(n)))
    beta = 0.0873

    #Start the central loop
    while not problem.final_target_hit and problem.evaluations < budget:
        
        #Perform recombination for a mu larger than 1
        if mu > 1:
            if rec == 'glob':
                if rec_type == 'discrete':
                    x_prime, sigma_prime, alpha_prime = global_discrete(x, sigma, alpha, lamb,n)
                elif rec_type == 'inter':
                    x_prime, sigma_prime, alpha_prime = global_inter(x, sigma, alpha, lamb)
            elif rec == 'pairs':
                if rec_type == 'discrete':
                    x_prime, sigma_prime, alpha_prime = pairs_discrete(x, sigma, alpha, lamb, n)
                elif rec_type == 'inter':
                    x_prime, sigma_prime, alpha_prime = pairs_inter(x, sigma, alpha, lamb)
        if mu == 1:
            x_prime, sigma_prime, alpha_prime = np.tile(x, (lamb, 1)), np.tile(sigma, (lamb, 1)), np.tile(alpha, (lamb, 1))


        #Mutate
        x_prime_prime, sigma_prime_prime, alpha_prime_prime = mutation_ind(x_prime, sigma_prime, alpha_prime, beta, tao_prime, tao)

        #Select 
        if sel == 'mu+lamb':
            #Stack the parents on the offspring
            x_tot = np.vstack((x, x_prime_prime))
            sigma_tot = np.vstack((sigma, sigma_prime_prime))
            alpha_tot = np.vstack((alpha, alpha_prime_prime))

        else:
            x_tot = x_prime_prime
            sigma_tot = sigma_prime_prime
            alpha_tot = alpha_prime_prime           

        #Evaluate the objectove function on all rows of the offspring matrix
        f = {i:problem(x_tot[i,:]) for i in range(len(x_tot))}
        f_sorted = sorted(f.items(), key=lambda item: item[1])
        indexes = [i[0] for i in f_sorted[:mu]]
        fopt = f_sorted[0][1]
        x = x_tot[indexes,:]
        sigma = sigma_tot[indexes,:]
        alpha = alpha_tot[indexes,:]
    return x, fopt



if __name__ == '__main__':

    ## Declarations of Ids, instances, and dimensions that the problems to be tested.
    problem_id = range(1,25)
    instance_id = range(1,26)
    dimension = [2, 5, 20]

    mu = 30
    lamb = 100

    '''
    The following parameters are already set on the best configuration we found. However. the user can change them if they want by using the
    suggested values in the comment.
    '''
        
    #Set it to "pairs" for pairwise recombination and "glob" for global recombination.
    rec = 'glob'
    #Set it to "inter" for intermediate and "discrete" for discrete
    rec_type = 'inter'
    #Set it to either "mu+lambda" or "mu,lambda"
    sel = 'mu,lamb'
    
    
    #Launch the algorithm
    logger = IOH_logger("./", "rot_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,"rot_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type,
                       "rot_"+sel+'_'+str(mu)+'_'+str(lamb)+'_'+rec+'_'+rec_type)

    for p_id in problem_id :
        for d in dimension :
            for i_id in instance_id:
                ## Getting the problem with corresponding id,dimension, and instance.
                f = IOH_function(p_id, d, i_id, suite="BBOB")
                #f.add_logger(logger)
                xopt, fopt = optimiz(f)
    logger.clear_logger()
