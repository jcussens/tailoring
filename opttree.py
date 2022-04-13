from gurobipy import *
from math import floor
import argparse
import numpy as np
import csv
import sys # for early exits in debugging
from itertools import combinations

def improve(posa,nega,rest,bvals,x,y,current_obj):
    """
    Given a current solution find the best zero value to change to -1 or 1
    """
    current_score = []
    for xi in x:
        cs = sum([xi[j] for j in posa]) - sum([xi[j] for j in nega])
        current_score.append(cs)
    best_obj = current_obj
    best_j = None
    best_sign = None
    best_b = None
    for j in rest:
        pos_scores = [cs + x[i][j] for i, cs in enumerate(current_score)]
        neg_scores = [cs - x[i][j] for i, cs in enumerate(current_score)]
        for scores, sign in (pos_scores,"pos"), (neg_scores,"neg"):
            for b in bvals:
                # get quality of policy with pos[j] newly set to 1 and threshold b
                obj = 0
                for i, yval in enumerate(y):
                    if scores[i] >= b+1:
                        obj += y[i]
                    if obj > best_obj:
                        best_j = j
                        best_b = b
                        best_obj = obj
                        best_sign = sign
    return best_j, best_sign, best_b, best_obj

def best_initial(y):
    treat_all_obj = sum(y)
    if treat_all_obj > 0:
        best_b = -1
        best_obj = treat_all_obj
    else:
        best_b = 0
        best_obj = 0
    return best_b, best_obj

def greedy_initial(x,y,p,bvals,maxnonzeros):
    posa = set()
    nega = set()
    rest = set(range(p))
    best_b, best_obj = best_initial(y)
    for steps in range(maxnonzeros):
        best_j, best_sign, best_b, best_obj = improve(posa,nega,rest,bvals,x,y,best_obj) 
        if best_j is None:
            break
        else:
            rest.remove(best_j)
            if best_sign == "pos":
                posa.add(best_j)
            else:
                nega.add(best_j)
    return posa, nega, best_b, best_obj
                    
def learn_tree(x,y,depth,x_names,solution_vector=None,report=False):

    p = len(x[0])
    n = len(x)

    model = Model()
    model.ModelSense = -1             # maximise objective (rather than the default of minimisation)
    #model.Params.CoverCuts = 2
    #model.Params.MIPFocus = 2
    
    # create variables for the hyperplane
    a = [] 
    for j in range(p):
        v = model.addVar(name="a#{0}".format(j),vtype=GRB.INTEGER,lb=-1,ub=1)
        a.append(v)

    # posa = []
    # nega = [] 
    # for j in range(p):
    #     v = model.addVar(name="posa#{0}".format(j),vtype=GRB.BINARY)
    #     posa.append(v)
    #     v = model.addVar(name="nega#{0}".format(j),vtype=GRB.BINARY)
    #     nega.append(v)

    # create variables indicating policy for each unit
    # note objective value is y[i]
    unit_policy = []
    for i in range(n):
        v = model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i])
        #v.BranchPriority = 10
        unit_policy.append(v)
        #unit_policy.append(model.addVar(name="up#{0}".format(i),vtype=GRB.CONTINUOUS,lb=0,ub=1,obj=y[i]))

    # create score variable for each unit
    #score = []
    maxvals = []
    for i in range(n):
        maxval = sum([x[i][j] for j in range(p)])
        #v = model.addVar(name="s#{0}".format(i),vtype=GRB.INTEGER,lb=-maxval,ub=maxval)
        #score.append(v)
        maxvals.append(maxval)

    # coefficient cannot be both -1 and +1 
    #for j in range(p):
    #    model.addConstr(posa[j] + nega[j] <= 1)

    maxnonzeros = 2
    
    #model.addConstr(quicksum(posa) + quicksum(nega) <= maxnonzeros)
    
    r = model.addVar(vtype=GRB.INTEGER,lb=0,ub=maxnonzeros)
    model.addGenConstrNorm(r, a, 2.0, "normconstr")

    
    # hyperplane threshold
    #b_bound = max(maxvals)
    b = model.addVar(name="b",vtype=GRB.INTEGER,lb=-maxnonzeros,ub=maxnonzeros-1)

    # score constraints
    #for i in range(n):
    #    model.addConstr(LinExpr([x[i][j] for j in range(p)],a) == score[i])

    # side of threshold determines unit policy
    for i in range(n):
        maxval = maxvals[i]
        # if up[i] = 0 then score[i] <= b
        #model.addConstr(LinExpr([x[i][j] for j in range(p)],posa) - LinExpr([x[i][j] for j in range(p)],nega) <= b + 2*maxval*unit_policy[i])
        #model.addConstr((unit_policy[i]==0) >> ((LinExpr([x[i][j] for j in range(p)],posa) - LinExpr([x[i][j] for j in range(p)],nega)) <= b))
        model.addConstr((unit_policy[i]==0) >> (LinExpr([x[i][j] for j in range(p)],a) <= b))
        # if up[i] = 1 then score[i] >= b+1
        #model.addConstr(LinExpr([x[i][j] for j in range(p)],posa) - LinExpr([x[i][j] for j in range(p)],nega) >= b + 1 - 2*maxval*(1-unit_policy[i]))
        model.addConstr((unit_policy[i]==1) >> (LinExpr([x[i][j] for j in range(p)],a) >= b+1))
        #model.addConstr((unit_policy[i]==1) >> ((LinExpr([x[i][j] for j in range(p)],posa) - LinExpr([x[i][j] for j in range(p)],nega)) >= b+1))


    # best_posa, best_nega, best_b, best_obj = greedy_initial(x,y,p,tuple(range(-maxnonzeros,maxnonzeros)),maxnonzeros)

    # model.NumStart = 1
    # model.Params.StartNumber = 0
    # b.Start = best_b
    # #for vs, best in (posa,best_posa), (nega,best_nega):
    # for vs, best in (a,best_posa), (a,best_nega):
    #     for j, v in enumerate(vs):
    #         if j in best:
    #             if best is best_posa
    #             v.Start = 1
    #         else:
    #             v.Start = 0
    # for i, xi in enumerate(x):
    #     score = 0
    #     for j in best_posa:
    #         score += xi[j]
    #     for j in best_nega:
    #         score -= xi[j]
    #     if score <= best_b:
    #         unit_policy[i].Start = 0
    #     else:
    #         unit_policy[i].Start = 1
            
    model.optimize()

    # posvs = []
    # negvs = []
    # for text, vs, lst in ('Positive',posa,posvs), ('Negative',nega,negvs):
    #     print('{0} covariates'.format(text))
    #     for j, v in enumerate(vs):
    #         if v.X > 0.5:
    #             lst.append(j)
    #             print(x_names[j])

    for j,v in enumerate(a):
        if v.X > 0.5 or v.X < -0.5:
            print(x_names[j],v.X)
    print('If score >= {0} then treat'.format(b.X+1))

    print('If all treated objective is {0}'.format(sum(y)))
    
    # for i in range(n):
    #     score = sum([x[i][j] for j in posvs]) - sum([x[i][j] for j in negvs])
    #     if unit_policy[i].X > 0.5:
    #         treat = True
    #         correct = (y[i]>=0)
    #     else:
    #         treat = False
    #         correct = (y[i]<0)
    #     print('{0}, score={1}, Treat={2}, yval={3}, correct={4}'.format(i,score,treat,y[i],correct)) 
        
        # for i in range(n):
        # if y[i] < 0:
        #     if unit_policy[i].X > 0.5:
        #         print('Treated {0} {1}'.format(i,y[i]))
        # else:
        #     if unit_policy[i].X < 0.5:
        #         print('Untreated {0} {1}'.format(i,y[i]))
    # Extract answer

    #For debugging
    #for v in model.getVars():
    #    print(v,v.X)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Learn an optimal policy tree from covariates and scores')
    parser.add_argument('data', help='file containing the CSV data (no commas in variables names allowed!)')
    parser.add_argument('--depth', default=2, type=int, help='(Maximum) depth of policy tree')
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
    
    subsample_size = 2000
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

    learn_tree(x,y,args.depth,x_names,solution_vector=None,report=True)
    
