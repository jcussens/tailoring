from gurobipy import *
from math import floor
import argparse
import numpy as np
import csv
import sys # for early exits in debuggin
from itertools import combinations

def get_best_j(x,y,p,indices):
    '''Find the best covariate to split on given a subset of the data
    '''
    scores = []
    for j in range(p):
        left = []
        right = []
        for i in indices:
            if x[i][j] == 0:
                left.append(i)
            else:
                right.append(i)
        score = max(0,sum([y[i] for i in left])) + max(0,sum([y[i] for i in right]))
        scores.append((score,j))
    scores.sort(reverse=True)
    bestj = scores[0][1]
    best_left = []
    best_right = []
    for i in indices:
        if x[i][bestj] == 0:
            best_left.append(i)
        else:
            best_right.append(i)
    return scores, bestj, best_left, best_right

def greedy(x,y,p,n,left_child,right_child,num_branch_nodes,num_nodes):
    '''Find a tree using greedy algorithm
    '''
    indices = [None]*num_nodes
    solution_vector = [None]*num_nodes
    indices[0] = list(range(n))
    for t in range(num_branch_nodes):
        _, bestj, best_left, best_right = get_best_j(x,y,p,indices[t])
        solution_vector[t] = bestj
        indices[left_child[t]] = best_left
        indices[right_child[t]] = best_right
    for t in range(num_branch_nodes,num_nodes):
        if sum([y[i] for i in indices[t]]) > 0:
            solution_vector[t] = 1
        else:
            solution_vector[t] = 0
    return solution_vector
            
def get_policy(xi,left_child,right_child,num_branch_nodes,solution_vector):
    t = 0
    while t < num_branch_nodes:
        if solution_vector[t] is None or xi[solution_vector[t]] == 0:
            t = left_child[t]
        else:
            t = right_child[t]
    return t, solution_vector[t]
        
def mycallback(model,where):
    '''Add lazy constraints given an integer solution
    '''
    #print(where)
    if where == GRB.Callback.MIPSOL:
        cutting = False
    elif (False and where == GRB.Callback.MIPNODE and
          model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
        cutting = True
    else:
        return
    a, node_policy, unit_policy, x, n, p, left_child, right_child, vals = model._info
    num_branch_nodes = len(a)
    #ncons = 0
    #print('Callback start')
    epsilon = 0.001
    for i in range(n):
        xi = x[i]
        valsi = vals[i]
        # what policy does the current tree assign to x[i] ?
        t = 0
        pathi = [0]
        pathvals = []
        cut = True
        while t < num_branch_nodes:
            if cutting:
                rightval = sum([model.cbGetNodeRel(a[t][j]) for j in valsi[1]])
            else:
                rightval = sum([model.cbGetSolution(a[t][j]) for j in valsi[1]])
            if rightval > 1 - epsilon:
                t = right_child[t]
                pathvals.append(valsi[1])
            elif rightval < epsilon:
                t = left_child[t]
                pathvals.append(valsi[0])
            else:
                cut = False
                break
            pathi.append(t)
        if not cut:
            # no cut from datapoint xi
            continue
        if cutting:
            tree_policyi = model.cbGetNodeRel(node_policy[t])
            sol_policyi = model.cbGetNodeRel(unit_policy[i])
        else:
            tree_policyi = model.cbGetSolution(node_policy[t])
            sol_policyi = model.cbGetSolution(unit_policy[i])
            # add a constraint if tree-determined policy for xi contradicts value in current solution
        if tree_policyi != sol_policyi:
            #print('Datapoint {0} should be {1} but is {2}'.format(i,tree_policyi,sol_policyi))
            #print([xi[j] for j in splits])
            #ncons += 1
            lexpr = LinExpr()
            #current_lhs = 0
            for k, t in enumerate(pathi[:-1]):
                lexpr.addTerms([-1]*len(pathvals[k]),[a[t][jj] for jj in pathvals[k]])
                #current_lhs -= model.cbGetSolution(a[t][splits[k]])
            #current_lhs -= model.cbGetSolution(node_policy[pathi[-1]])
            if tree_policyi > sol_policyi:
                lexpr.add(node_policy[pathi[-1]],-1)
                lexpr.add(unit_policy[i],1)
                #current_lhs += model.cbGetSolution(unit_policy[i])
            else:
                lexpr.add(node_policy[pathi[-1]],1)
                lexpr.add(unit_policy[i],-1)
                #current_lhs -= model.cbGetSolution(unit_policy[i])
            #print('new cons:', lexpr, '>=', 1-len(pathi))
            if cutting:
                model.cbCut(lexpr, GRB.GREATER_EQUAL, 1-len(pathi))
            else:
                model.cbLazy(lexpr, GRB.GREATER_EQUAL, 1-len(pathi))
    #print('Callback done')
    #print('New incumbent generated {0} constraints'.format(ncons))
                    
def learn_tree(x,y,depth,x_names,solution_vector=None,report=False):

    p = len(x[0])
    n = len(x)

    #print(n,p)
    
    num_nodes = 2**(depth+1)-1
    num_branch_nodes = floor(num_nodes/2)
    num_leaf_nodes = num_nodes - num_branch_nodes

    # Construct left-branch and right-branch ancestors for each node
    # Use breadth-first node numbering
    # Root is node 0
    # Next level has nodes 1 and 2, left and right child of 0 resp.
    # Next level has nodes 3,4 (children of node 1) and 5,6 (children of node 2)
    # and so on

    level = 0
    frontier = [0]
    current_node = 0
    al = {0:set()}
    ar = {0:set()}
    parent = {}
    left_child = {}
    right_child = {}
    at_depth = [0]
    while level < depth:
        new_frontier = []
        for node in frontier:
            left = current_node+1
            right = current_node+2
            at_depth.extend([level+1,level+1])
            parent[left] = node
            parent[right] = node
            left_child[node] = left
            right_child[node] = right
            new_frontier.extend([left,right])
            al[left] = al[node] | set([node])
            ar[left] = ar[node]
            al[right] = al[node]
            ar[right] = ar[node] | set([node])
            current_node = right
        frontier = new_frontier
        level += 1

    # Debugging/check
    #print('left',al)
    #print('right',ar)

    model = Model()
    model.ModelSense = -1             # maximise objective (rather than the default of minimisation)
    model.Params.LazyConstraints = 1  # constraints added lazily
    model.Params.PreCrush = 1         # constraints added lazily
    #model.Params.BranchDir = 1
    #model.Params.MIPFocus = 3
    model.Params.PreSOS1Encoding = 3
    
    # create variables for the tree structure
    # a[t][j]=1 iff we split on j at node t. NB a[j][t] in paper
    # if x[i][j]=1 and datapoint i reaches node t and a[t][j]=1 then it goes right
    # if x[i][j]=0 and datapoint i reaches node t and a[t][j]=1 then it goes left
    a = [] 

    for t in range(num_branch_nodes):
        a.append([])
        for j in range(p):
            v = model.addVar(name="a#{0}#{1}".format(t,j),vtype=GRB.BINARY)
            v.BranchPriority = 100-at_depth[t]
            a[t].append(v)
            
    # create variables indicating policy for each leaf
    node_policy = [None]*num_branch_nodes
    for t in range(num_branch_nodes,num_nodes):
        node_policy.append(model.addVar(name="np#{0}".format(t),vtype=GRB.BINARY))
        #node_policy.append(model.addVar(name="np#{0}".format(t),vtype=GRB.CONTINUOUS,lb=0,ub=1))

    # create variables indicating policy for each unit
    # note objective value is y[i]
    unit_policy = []
    for i in range(n):
        v = model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i])
        #v.BranchPriority = 10
        unit_policy.append(v)
        #unit_policy.append(model.addVar(name="up#{0}".format(i),vtype=GRB.CONTINUOUS,lb=0,ub=1,obj=y[i]))

    # constraints for the tree structure

    scores, bestj, best_left, best_right = get_best_j(x,y,p,list(range(n)))

    for (_,j) in scores[100:]:
        for t in range(num_branch_nodes):
            model.addConstr(a[t][j] == 0)

    weights = [s[0] for s in scores]
    for t in range(num_branch_nodes):
        # Split on exactly one covariate in each node
        model.addConstr(quicksum(a[t]) == 1,name="a{0}".format(t))
        aa = [a[t][s[1]] for s in scores]
        model.addSOS(GRB.SOS_TYPE1, aa, weights)


    # force neighbouring leaves to have opposite policies
    # for t in range(num_branch_nodes,num_nodes):
    #     # if t is odd ...
    #     if t % 2 == 1:
    #         model.addConstr(node_policy[t] + node_policy[t+1] == 1,name="py{0}".format(t))

    # for each datapoint find covariates where it = 0 and those where it = 1
    vals = []
    for i in range(n):
        val_zero = []
        val_one = []
        for j, val in enumerate(x[i]):
            if val == 0:
                val_zero.append(j)
            else:
                val_one.append(j)
        vals.append((frozenset(val_zero),frozenset(val_one)))

    # for i1,i2 in combinations(range(n),2):
    #     if x[i1] == x[i2]:
    #         print('Datapoints {0} and {1} cannot be distinguished'.format(i1,i2))
    #         model.addConstr(
    # sys.exit()
        
    # For each pair of datapoints with opposite policies, state that they get
    # the same policy (due to being in the same leaf) if no covariate distinguishing them
    # is in the tree
    neg_indices = []
    pos_indices = []
    for i in range(n):
        if y[i] >= 0:
            pos_indices.append(i)
            #model.addConstr(unit_policy[i] == 1)
        else:
            #model.addConstr(unit_policy[i] == 0)
            neg_indices.append(i)

    lazyval = 3
            
    diffconslim = 5
    for negi in neg_indices:
        zeroval_covariate_negi = vals[negi][0]
        for posi in pos_indices:
            diff = zeroval_covariate_negi ^ vals[posi][0]
            #print('To distinguish datapoint {0} and {1} one of these must be in tree {2}'.format(negi,posi,diff))
            if len(diff) <= diffconslim:
                lexpr = LinExpr()
                for t in range(num_branch_nodes):
                    lexpr.addTerms([1]*len(diff),[a[t][j] for j in diff])
                cons = model.addConstr(lexpr + unit_policy[negi] - unit_policy[posi] >= 0,name="dg1{0}{1}".format(negi,posi))
                cons.Lazy = lazyval 
                cons = model.addConstr(lexpr - unit_policy[negi] + unit_policy[posi] >= 0,name="dg1{0}{1}".format(negi,posi))
                cons.Lazy = lazyval

    samediffconslim = 2
    for indices in neg_indices, pos_indices:
        for i, idx1 in enumerate(indices):
            zeroval_covariate_idx1 = vals[idx1][0]
            for idx2 in indices[i+1:]:
                diff = zeroval_covariate_idx1 ^ vals[idx2][0]
                if len(diff) <= samediffconslim:
                    lexpr = LinExpr()
                    for t in range(num_branch_nodes):
                        lexpr.addTerms([1]*len(diff),[a[t][j] for j in diff])
                    cons = model.addConstr(lexpr + unit_policy[idx1] - unit_policy[idx2] >= 0)
                    cons.Lazy = lazyval
                    cons = model.addConstr(lexpr - unit_policy[idx1] + unit_policy[idx2] >= 0)
                    cons.Lazy = lazyval
            
            
    # add information as attributes to model, so available in callback
    model._info = a, node_policy, unit_policy, x, n, p, left_child, right_child, vals

    # get greedy solution
    solution_vector = greedy(x,y,p,n,left_child,right_child,num_branch_nodes,num_nodes)
    
    if solution_vector is not None:
        model.NumStart = 1
        model.Params.StartNumber = 0
        for t in range(num_branch_nodes):
            for j in range(p):
                a[t][j].Start = 0
            a[t][solution_vector[t]].Start = 1
        for t in range(num_branch_nodes,num_nodes):
            node_policy[t].Start = solution_vector[t]
        for i in range(n):
            _, unit_policy[i].Start = get_policy(x[i],left_child,right_child,num_branch_nodes,solution_vector)
    
    model.optimize(mycallback)

    # Extract answer

    solution_vector = []
    
    # The tree

    if report:
        branches = []
        for (node,node_parent) in parent.items():
            branches.append((node_parent,node))
        branches.sort()
        print('Tree structure:')
        for node_parent,node in branches:
            print('{0}->{1}'.format(node_parent,node))
        print()

    for t in range(num_branch_nodes):
        for j in range(p):
            if a[t][j].X > 0.5:
                solution_vector.append(j)
                if report:
                    print("Node {0} splits on covariate {1}".format(t,x_names[j]))
                break
    if report:
        print()

    for t in range(num_branch_nodes,num_nodes):
        if node_policy[t].X > 0.5:
            np = 1
        else:
            np = 0
        solution_vector.append(np)
        if report:
            print("Leaf node {0} chooses policy {1}".format(t,np))
    if report:
        print()

    if False and report:
        for i in range(n):
            t, pi = get_policy(x[i],left_child,right_child,num_branch_nodes,solution_vector)
            print("Unit {0} (obj={3}) ends up in leaf node {1} and is assigned policy {2}".format(i,t,pi,y[i]))
            
    return solution_vector
        
    #for i in range(n):
    #    for t in range(num_branch_nodes,num_nodes):
    #        if z[i][t].X > 0.5:
    #            print("Unit {0} (obj={3}) ends up in leaf node {1} and is assigned policy {2}".format(i,t,unit_policy[i].X,y[i]))
    #            break

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
    
    subsample_size = 200
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
    
