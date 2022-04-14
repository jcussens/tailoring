from gurobipy import *
from math import floor
import argparse
import numpy as np
import csv
import sys # for early exits in debugging
from itertools import combinations
                    
def learn_tree(x,y,x_names,report=False):

    p = len(x[0])
    n = len(x)

    model = Model()
    model.ModelSense = -1             # maximise objective (rather than the default of minimisation)
    #model.Params.CoverCuts = 2
    #model.Params.MIPFocus = 1
    
    # create variables for the hyperplane
    a = [] 
    for j in range(p):
        v = model.addVar(name="a#{0}".format(j),vtype=GRB.CONTINUOUS,lb=-1,ub=1)
        a.append(v)

    # create variables indicating policy for each unit
    # note objective value is y[i]
    unit_policy = []
    for i in range(n):
        v = model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i])
        #v.BranchPriority = 10
        unit_policy.append(v)

    maxnonzeros = 5
    
    r = model.addVar(vtype=GRB.INTEGER,lb=0,ub=maxnonzeros)
    model.addGenConstrNorm(r, a, 0.0, "normconstr")

    b = 0
    
    # side of threshold determines unit policy
    for i in range(n):
        if y[i] < 0:
            model.addConstr(LinExpr([x[i][j] for j in range(p)],a) <= b + (maxnonzeros-b)*unit_policy[i])
        else:
            model.addConstr(LinExpr([x[i][j] for j in range(p)],a) >= b + 1 - (maxnonzeros+b+1)*(1-unit_policy[i]))
            
    model.optimize()

    for j,v in enumerate(a):
        if v.X > 0.001 or v.X < -0.001:
            print(x_names[j],v.X)
    print('If score >= {0} then treat'.format(b))

    print('If all treated objective is {0}'.format(sum(y)))
    

    #For debugging
    #for v in model.getVars():
    #    print(v,v.X)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Learn an optimal policy tree from covariates and scores')
    parser.add_argument('data', help='file containing the CSV data (no commas in variables names allowed!)')
    parser.add_argument('--exclude', default="A,Y", help='Comma separated string of covariates to exclude')
    parser.add_argument('--scores', default="scores.DML,scores.DR,scores.CF", help='Comma separated string of fields to interpret as scores')
    parser.add_argument('--score', default="scores.DML", help='Score to use to build tree')
    
    args = parser.parse_args()

    data = np.genfromtxt(args.data,delimiter=',',skip_header=True)
    header = open(args.data).readline().rstrip()

    #Use CSV reader to correctly parse header line
    for line in csv.reader([header],skipinitialspace=True):
        names = line
        break

    # construct x and y
    
    x_names_to_exclude = frozenset(args.exclude.split(",") + args.scores.split(","))
    x_indices_to_exclude = []
    x_names = []
    for i, name in enumerate(names):
        if name in x_names_to_exclude:
            x_indices_to_exclude.append(i)
        else:
            x_names.append(name)
        if name == args.score:
            y_index = i
    x = np.delete(data,x_indices_to_exclude,axis=1)
    x = x.astype(np.uint32)
    y = data[:,y_index]


    # choose a subsample
    
    subsample_size = 200000
    x = x[:subsample_size]
    y = y[:subsample_size]


    # collapse identical covariate vectors
    
    x, indices  = np.unique(x,axis=0,return_inverse=True)
    newy = np.zeros(len(x))
    for i, idx in enumerate(indices):
        newy[idx] += y[i]
    y = newy
    
    # switch to lists
    
    x = x.tolist()
    y = y.tolist()

    learn_tree(x,y,x_names)
    
