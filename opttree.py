from gurobipy import *
from math import floor

def ladd(xx,yy):
    '''add two lists together
    '''
    return [x + yy[i] for i, x in enumerate(xx)]

depth = 3

#fake covariates
x = [
    [0.2,0.1,0.0,0.3],
    [0.3,0.5,0.2,0.4],
    [0.1,0.2,0.3,0.4],
    [0.5,0.6,0.7,0.8],
    [0.9,1.0,0.8,0.7],
    [0.3,0.3,0.3,0.3],
    [0.2,0.2,0.0,0.9],
    [0.15,0.1,0.3,0.3],
    [0.7,0.9,0.2,0.3],
    [0.2,0.3,0.4,0.5]]

# fake estimates of benefit of treatment
y = [-2,3,6,0.2,-12,-3,-5,-6,7,5]

x = x*300
y = y*300

p = len(x[0])
n = len(x)


# find each "epsilon_j"
epsilon = []
for j in range(p):
    xj = [x[i][j] for i in range(n)]
    xj.sort()
    epsilon.append(min([xj[i+1]-xj[i] for i in range(n-1) if xj[i+1]-xj[i] > 0]))
e_max = max(epsilon)

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
while level < depth:
    new_frontier = []
    for node in frontier:
        left = current_node+1
        right = current_node+2
        parent[left] = node
        parent[right] = node
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

# variables for the tree structure

a = [] # a[t][j]=1 iff we split on j at node t. NB a[j][t] in paper
d = [] # d[t]=1 iff there is a split at node t
b = [] # b[t] is the splitting threshold at node t

for t in range(num_branch_nodes):
    d.append(model.addVar(name="d{0}".format(t),vtype=GRB.BINARY))
    b.append(model.addVar(name="b{0}".format(t),vtype=GRB.CONTINUOUS,lb=0,ub=1))
    a.append([])
    for j in range(p):
        a[t].append(model.addVar(
            name="a{0}{1}".format(t,j),vtype=GRB.BINARY))

# constraints for the tree structure

for t in range(num_branch_nodes):
    
    # Constraint (2) in the paper
    model.addConstr(quicksum(a[t]) == d[t])

    #Currently requiring a split at each node
    model.addConstr(d[t] == 1)

    # Constraint (3) in the paper
    model.addConstr(b[t] <= d[t])

for t in range(1,num_branch_nodes):

    # Constraint (5) in the paper
    model.addConstr(d[t] <= d[parent[t]])

# variables for datapoints
    
z = [] # z[i][t] =  1 iff datapoint i is "in" leaf node t
for i in range(n):
    z.append([None]*num_branch_nodes)
    for t in range(num_branch_nodes,num_nodes):
        z[i].append(model.addVar(name="z{0}{1}".format(i,t),vtype=GRB.BINARY))

# constraints assigning datapoints to nodes

for i in range(n):
    # each datapoint in exactly one leaf
    model.addConstr(quicksum(z[i][num_branch_nodes:]) == 1)

# finally the splits!
# NB there is a typo in the paper in (24) where they have b_t when it should be b_m
# Bertsimas and Dunn paper is wrong! (11) is correct, but (12) is wrong when there is no split and where everything should go left
for i in range(n):
    for t in range(num_branch_nodes,num_nodes):
        for m in ar[t]:
            model.addConstr(LinExpr(x[i],a[m]), GRB.GREATER_EQUAL, b[m] - (1-z[i][t]))
        for m in al[t]:
            model.addConstr(LinExpr(ladd(x[i],epsilon),a[m]), GRB.LESS_EQUAL, b[m] + (1+e_max)*(1-z[i][t]))

# NOT FROM THE PAPER

# create variables indicating policy for each leaf
node_policy = [None]*num_branch_nodes
for t in range(num_branch_nodes,num_nodes):
    node_policy.append(model.addVar(name="np{0}".format(t),vtype=GRB.BINARY))

# create variables indicating policy for each unit
# note objective value is y[i]
unit_policy = []
for i in range(n):
    unit_policy.append(model.addVar(name="up{0}".format(i),vtype=GRB.BINARY,obj=y[i]))

#for i in range(n):
#    model.addConstr(unit_policy[i] == quicksum([z[i][t] * node_policy[t] for t in range(num_branch_nodes,num_nodes)]))

for i in range(n):
    for t in range(num_branch_nodes,num_nodes):
        model.addConstr(node_policy[t] + (1-z[i][t]) + (1-unit_policy[i]) >= 1)
        model.addConstr((1-node_policy[t]) + (1-z[i][t]) + unit_policy[i] >= 1)
    
# optimisation
# neighbouring leaves have opposite policies
for t in range(num_branch_nodes,num_nodes):
    # if t is odd ...
    if t % 2 == 1:
        model.addConstr(node_policy[t] + node_policy[t+1] == 1)
# optional resource constraint, at most 2 allowed to get treatment!
#model.addConstr(quicksum(unit_policy) <= 2)
    
# last constraint set contains non-convex quadratic constraint, so have to tell Gurobi we're OK with that
# they are bilinear constraints so should not cause too much trouble.
#model.params.NonConvex = 2

model.ModelSense = -1      # maximise objective (rather than the default of minimisation)

model.optimize()

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
    if d[t].X < 0.5:
        print('No split at node {0}'.format(t))
    else:
        for j in range(p):
            if a[t][j].X > 0.5:
                break
        print("Node {0} splits on covariate {1} at threshold of {2}".format(t,j,b[t].X))
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
