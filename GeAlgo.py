'''This is the code for genetic algorithm. 3 Parts are included
- the main section, the selection part, and the genchild part
'''
import numpy as np
import numpy.matlib as npm
from scipy.special import comb
import csv

def Selection(e, g, pop):
    I = np.argsort(e)
    etemp = e[I]
    gtemp = g[I,:]

    e0 = np.zeros(pop)
    g0 = np.zeros((pop,len(g[1,:])))
    e0[0] = etemp[0]
    g0[0,:] = gtemp[0,:]

    p = 0.3
    psbly = np.empty(len(e) - 1)
    for i in range(1, len(e)):
        psbly[i - 1] = p * (1 - p) ** i

    psbly = psbly / sum(psbly)
    psblymulti = 0
    psblytrue = np.zeros(len(psbly) + 1)
    for i in range(len(psbly)):
        psblymulti = psblymulti + psbly[i]
        psblytrue[i + 1] = psblymulti

    i = 1
    list = []
    while i < pop:
        Randprocess = np.random.rand()
        for j in range(1, len(e)):
            if (Randprocess <= psblytrue[j]) and (Randprocess >= psblytrue[j - 1]):
                if j in list:
                    i = i - 1
                    break
                else:
                    # disp('New list!');
                    e0[i] = etemp[j]
                    g0[i, :] = gtemp[j, :]
                    list.append(j)

        i = i + 1

    return e0, g0

def GenChild(parents, lb, ub, mut):

    bs = ub - lb
    count = 0
    s = parents.shape

    children = np.zeros([int(comb(s[0], 2)), s[1]])

    # Crossover
    for i in range(s[0] - 1):
        for j in range(i + 1, s[0]):
            count += 1
            for k in range(s[1]):
                rnd = np.random.rand()
                if rnd < 1/3:
                    children[count - 1, k] = parents[i, k]
                elif rnd < 2/3:
                    children[count - 1, k] = parents[j, k]
                else:
                    children[count - 1, k] = (parents[i, k ] + parents[j, k]) / 2

    # Mutation
    children = children + np.random.normal(0, mut, size = (count, s[1])) * npm.repmat(bs, count, 1)

    # Check and correct boundary
    for i in range(count):
        for j in range(s[1]):
            if children[i, j] > ub[j]:
                children[i, j] = ub[j]
            if children[i, j] < lb[j]:
                children[i, j] = lb[j]


    return children

def GeAlgo( fun, lb, ub, algorithm_parameters):
    pop = algorithm_parameters['population_size']
    it_tot = algorithm_parameters['max_num_iteration']
    mut = algorithm_parameters['mutation_probability']

    lparameter = len(lb)

    count = int(comb(pop, 2))

    g0 = np.zeros([pop, lparameter])
    split_l = int(np.floor(pop / 2))
    split_u = split_l + 1
    num_l = int(np.floor(pop / 2))
    num_u = int(np.floor(pop / 2 + 0.5 ))
    half_spread = (ub - lb) / 2

    g0[0 : split_l, :] = np.random.rand(num_l, lparameter) * npm.repmat(half_spread, num_l, 1) + npm.repmat(lb, num_l, 1)
    g0[split_u - 1: ,:] = np.random.rand(num_u, lparameter) * npm.repmat(half_spread, num_u, 1) + npm.repmat(lb + half_spread, num_u, 1)

    record_file = "record_gealg.csv"
    record_title = ['Generation', 'a', 'b', 'af', 'bf', 'as', 'bs', 'e']
    with open(record_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(record_title)

    for j in range(lparameter):
        g0[:, j] = g0[np.random.permutation(pop), j]

    print('*** g0 Generated ')

    e0 = np.zeros(pop)
    for k in range(pop):
        e0[k] = fun(g0[k, :])

    print('*** g0 Generation Finished')

    beste = np.zeros(it_tot)
    gen_term = algorithm_parameters['mutation_change_generation']
    tol_term = 0.01
    mut_factor = algorithm_parameters['mutation_change_factor']
    gen_count = 0

    for i in range(it_tot):
        print('*** Start of Generation # ' + str(i))
        gen_count = gen_count + 1

        g1 = GenChild(g0, lb, ub, mut)

        e1 = np.zeros(count)
        for k in range(count):
            e1[k] = fun(g1[k, :])

        g = np.concatenate((g1, g0))
        e = np.concatenate((e1, e0))
        e0, g0 = Selection(e, g, pop)

        beste[i] = e0[0]

        record_result = np.concatenate((g0, np.expand_dims(e0, axis = 1)), axis = 1)
        record_index  = npm.repmat(np.array(i), pop, 1)
        record_result = np.concatenate((record_index, record_result), axis = 1)
        with open(record_file, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerows(record_result)

        if i >= gen_term:
            if (np.abs(beste[i] - beste[i - gen_term]) < tol_term) and (gen_count >= gen_term):
                print('Early termination at generation ' + str(i))
                mut = mut / mut_factor
                print('Mutation is changed to ' + str(mut))
                gen_count = 0


    result = {'para': g0[0, :],
              'obj':  e0[0]
              }

    return result

# def fun(x):
#     return (x[0]-x[1]) ** 2 + (x[2]-x[3]) ** 2 + (x[4]-x[5]) ** 2
#
# lb = np.array([2,1,3,1,3,1])
# ub = np.array([5,50,6,50,6,50])
#
# algorithm_parameter = {'max_num_iteration': 30,
#                    'population_size':4,
#                    'mutation_probability':0.1,
#                     'mutation_change_generation': 5,
#                     'mutation_change_factor': 1.3}
#
# r = GeAlgo(fun, lb, ub, algorithm_parameter)
# print(r['para'], r['obj'])



