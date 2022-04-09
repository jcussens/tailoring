from gurobipy import *
from math import floor
import argparse
import numpy as np
import csv


def mycallback(model,where):
    '''Add lazy constraints given an integer solution
    '''
    #print(where)
    if where == GRB.Callback.MIPSOL:
        a, node_policy, unit_policy, x, n, p, left_child, right_child = model._info
        num_branch_nodes = len(a)
        ncons = 0
        #print('Callback start')
        for i in range(n):
            xi = x[i]
            # what policy does the current tree assign to x[i] ?
            t = 0
            pathi = [0]
            splits = []
            while t < num_branch_nodes:
                for j in range(p):
                    if model.cbGetSolution(a[t][j]) > 0.5:
                        splits.append(j)
                        break
                if xi[j] == 0:
                    t = left_child[t]
                else:
                    t = right_child[t]
                pathi.append(t)
            if model.cbGetSolution(node_policy[t]) > 0.5:
                tree_policyi = 1
            else:
                tree_policyi = 0

            if model.cbGetSolution(unit_policy[i]) > 0.5:
                sol_policyi = 1
            else:
                sol_policyi = 0
                
            # add a constraint if tree-determined policy for xi contradicts value in current solution
            if tree_policyi != sol_policyi:
                #print('Datapoint {0} should be {1} but is {2}'.format(i,tree_policyi,sol_policyi))
                #print([xi[j] for j in splits])
                ncons += 1
                lexpr = LinExpr()
                #current_lhs = 0
                for k, t in enumerate(pathi[:-1]):
                    lexpr.add(a[t][splits[k]],-1)
                    #current_lhs -= model.cbGetSolution(a[t][splits[k]])
                #current_lhs -= model.cbGetSolution(node_policy[pathi[-1]])
                if tree_policyi == 1:
                    lexpr.add(node_policy[pathi[-1]],-1)
                    lexpr.add(unit_policy[i],1)
                    #current_lhs += model.cbGetSolution(unit_policy[i])
                else:
                    lexpr.add(node_policy[pathi[-1]],1)
                    lexpr.add(unit_policy[i],-1)
                    #current_lhs -= model.cbGetSolution(unit_policy[i])
                #print('new cons:', lexpr, '>=', 1-len(pathi))
                model.cbLazy(lexpr, GRB.GREATER_EQUAL, 1-len(pathi))
        #print('Callback done')
        #print('New incumbent generated {0} constraints'.format(ncons))
                    
def learn_tree(x,y,depth):

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
    while level < depth:
        new_frontier = []
        for node in frontier:
            left = current_node+1
            right = current_node+2
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

    
    # create variables for the tree structure
    # a[t][j]=1 iff we split on j at node t. NB a[j][t] in paper
    # if x[i][j]=1 and datapoint i reaches node t and a[t][j]=1 then it goes right
    # if x[i][j]=0 and datapoint i reaches node t and a[t][j]=1 then it goes left
    a = [] 

    for t in range(num_branch_nodes):
        a.append([])
        for j in range(p):
            a[t].append(model.addVar(
                name="a#{0}#{1}".format(t,j),vtype=GRB.BINARY))

    # create variables indicating policy for each leaf
    node_policy = [None]*num_branch_nodes
    for t in range(num_branch_nodes,num_nodes):
        node_policy.append(model.addVar(name="np#{0}".format(t),vtype=GRB.BINARY))

    # create variables indicating policy for each unit
    # note objective value is y[i]
    unit_policy = []
    for i in range(n):
        unit_policy.append(model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i]))

    # constraints for the tree structure

    for t in range(num_branch_nodes):
        # Split on exactly one covariate in each node
        model.addConstr(quicksum(a[t]) == 1)

    # force neighbouring leaves to have opposite policies
    for t in range(num_branch_nodes,num_nodes):
        # if t is odd ...
        if t % 2 == 1:
            model.addConstr(node_policy[t] + node_policy[t+1] == 1)


    # add information as attributes to model, so available in callback
    model._info = a, node_policy, unit_policy, x, n, p, left_child, right_child
        
    model.optimize(mycallback)

    # Extract answer

    # The tree

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
                print("Node {0} splits on covariate {1}".format(t,j))
                break
    print()

    for t in range(num_branch_nodes,num_nodes):
        print("Leaf node {0} chooses policy {1}".format(t,node_policy[t].X))
    print()

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

    x_names_to_exclude = frozenset(args.exclude.split(",") + args.scores.split(","))
    x_indices_to_exclude = []
    for i, name in enumerate(names):
        if name in x_names_to_exclude:
            x_indices_to_exclude.append(i)
        if name == args.score:
            y_index = i

    x = np.delete(data,x_indices_to_exclude,axis=1)
    y = data[:,y_index]

    #x, indices, counts = np.unique(x,axis=0,return_index=True,return_counts=True)

    # switch to lists
    
    x = x.tolist()[:100]
    y = y.tolist()[:100]

    # newy = []
    # for i in indices:
    #     newy.append(y[i]*counts[i])
    # y = newy
        
    learn_tree(x,y,args.depth)
    
    #print(x)

    #print(y)

