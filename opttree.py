from gurobipy import *
import argparse
import numpy as np
import csv

class Node:

    def __init__(self,covariate=None,name=None,left=None,right=None,mipvars=None,
                 mother=None,indices=None,depth=None,policy=None):

        self._covariate = covariate
        self._name = name
        self._left = left
        self._right = right
        self._mipvars = mipvars
        self._mother = mother
        self._indices = indices
        self._depth = depth
        self._policy = policy

    def treestr(self):
        res = str(self)
        if self._left is not None:
            res += self._left.treestr()
        if self._right is not None:
            res += self._right.treestr()
        return res

    def __str__(self):
        res = '{0}\n'.format(self._name)
        if self._covariate is None:
            res += 'policy   : {0}\n'. format(self._policy)
        else:
            res += 'covariate: {0}\n'. format(self._covariate)
            res += 'left     : {0}\n'. format(self._left._name)
            res += 'right    : {0}\n'. format(self._right._name)
        return res + '\n'

    def get_path_from_root(self):
        mother = self._mother
        if mother is None:
            return []
        else:
            nodesdecisions = mother.get_path_from_root()
            if mother._left is self:
                decision = 0
            else:
                decision = 1
        return nodesdecisions + [(mother,decision)]


    def apply_solution(self):
        if self._left is not None:
            # if a left child then should always be a right child
            assert self._right is not None
            for j, v in enumerate(self._mipvars):
                if v.X > 0.5:
                    self._covariate = j
                    break
            self._left.apply_solution()
            self._right.apply_solution()
        else:
            if self._mipvars[0].X > 0.5:
                self._policy = 1
            else:
                self._policy = 0
            
    def addmipvars(self,model,p,x_names):
        if self._left is not None:
            # if a left child then should always be a right child
            assert self._right is not None
            self._mipvars = []
            for j in range(p):
                v = model.addVar(name="a#{0}#{1}".format(self._name,x_names[j]),vtype=GRB.BINARY)
                v.BranchPriority = 10-self._depth
                self._mipvars.append(v)
            self._left.addmipvars(model,p,x_names)
            self._right.addmipvars(model,p,x_names)
        else:
            v = model.addVar(name="np#{0}".format(self._name),vtype=GRB.BINARY)
            v.BranchPriority = 100
            self._mipvars = [v]
    
    def addmipcons(self,model,x,upvars,lazy,siblingcons,leftzerocons,usesos):
        if self._left is not None:
            # if a left child then should always be a right child
            assert self._right is not None
            model.addConstr(quicksum(self._mipvars) == 1)
            if usesos:
                model.addSOS(GRB.SOS_TYPE1, self._mipvars)
            self._left.addmipcons(model,x,upvars,lazy,siblingcons,leftzerocons,usesos)
            self._right.addmipcons(model,x,upvars,lazy,siblingcons,leftzerocons,usesos)
        else:
            npvar = self._mipvars[0]
            if self._policy is not None:
                model.addConstr(npvar == self._policy)


            if self._mother._left is self:
                if leftzerocons:
                    # just fix policy to 0 and set sibling's policy to 1
                    model.addConstr(npvar == 0)
                    model.addConstr(self._mother._right._mipvars[0] == 1)
                elif siblingcons:
                    # if a left child state that it has opposite policy from sibling
                    model.addConstr(npvar + self._mother._right._mipvars[0] == 1)
                
            nodesdecisions = self.get_path_from_root()
            for i, xi in enumerate(x):
                lexpr = LinExpr()
                for node, decision in nodesdecisions:
                    avars = [v for j, v in enumerate(node._mipvars) if xi[j] == decision]
                    lexpr.addConstant(1)
                    lexpr.addTerms([-1]*len(avars),avars)
                if self._policy is None:
                    cons = model.addConstr(lexpr + npvar - upvars[i] >= 0)
                    cons.Lazy = lazy
                    cons = model.addConstr(lexpr - npvar + upvars[i] >= 0)
                    cons.Lazy = lazy
                elif self._policy == 0:
                    cons = model.addConstr(lexpr + npvar - upvars[i] >= 0)
                    cons.Lazy = lazy
                else:
                    cons = model.addConstr(lexpr - npvar + upvars[i] >= 0)
                    cons.Lazy = lazy
                    

    def maketree(self,depth):
        if depth == 0:
            return
        left_node = Node(mother=self,depth=self._depth+1)
        right_node = Node(mother=self,depth=self._depth+1)
        self.set_children(left_node,right_node)
        left_node.maketree(depth-1)
        right_node.maketree(depth-1)
            
    def set_covariate(self,covariate):
        self._covariate = covariate

    def set_children(self,left_child,right_child):
        self._left = left_child
        self._right = right_child

    def set_policy(self,policy):
        '''Make a node a leaf node with given policy
        '''
        self._covariate = None
        self._left = None
        self._right = None
        self._policy = policy

    def score(self):
        if self._policy == 0:
            return 0
        elif self._policy == 1:
            return sum([y[i] for i in self._indices])
        else:
            return self._left.score() + self._right.score()
            
        
    def depth_first_enumerate(self,idx,x_names):
        self._name = str(idx)
        if self._covariate is not None:
            self._covariate = x_names[self._covariate]
        next = idx + 1
        if self._left is not None:
            next = self._left.depth_first_enumerate(next,x_names)
        if self._right is not None:
            next = self._right.depth_first_enumerate(next,x_names)
        return next

    def split_on_j(self,splits,j):
        left_indices =  splits[j][0] & self._indices
        right_indices = splits[j][1] & self._indices
        return Node(mother=self,indices=left_indices,depth=self._depth+1), Node(mother=self,indices=right_indices,depth=self._depth+1)
        
    def find_best_policy_and_score(self,y):
        score_for_1 = sum([y[i] for i in self._indices]) 
        if score_for_1 > 0:
            return 1, score_for_1
        else:
            return 0, 0

    def tree_greedy(self,x,splits,y,treedepth):
        '''
        Find a policy tree for the given subset of the data and depth
        using the greedy algorithm
        '''
        if self._depth == treedepth:
            best_policy, _ = self.find_best_policy_and_score(y)
            self.set_policy(best_policy)
        else:
            bestj, left_node, right_node = self.find_best_covariate(x,splits,y)
            self.set_covariate(bestj)
            self.set_children(left_node,right_node)
            left_node.tree_greedy(x,splits,y,treedepth)
            right_node.tree_greedy(x,splits,y,treedepth)


    def tree_optimal(self,x,splits,y,treedepth):
        '''
        Find an optimal policy tree for the given subset of the data and depth
        '''
        if self._depth == treedepth:
            best_policy, best_score = self.find_best_policy_and_score(y)
            self.set_policy(best_policy)
        else:
            p = len(x[0])  # assume non-empty data
            bestj = None
            for j in range(p):
                left_node, right_node = self.split_on_j(splits,j)
                best_left_score = left_node.tree_optimal(x,splits,y,treedepth)
                best_right_score = right_node.tree_optimal(x,splits,y,treedepth)
                score = best_left_score + best_right_score
                if bestj is None or score > best_score:
                    bestj = j
                    best_score = score
                    best_left_node = left_node
                    best_right_node = right_node
            self.set_covariate(bestj)
            self.set_children(best_left_node,best_right_node)
        return best_score

        
    def find_best_covariate(self,x,splits,y):
        '''Find the best covariate to split on
        Return it together with the two nodes it creates
        '''
        bestj = None
        p = len(x[0])  # assume non-empty data
        for j in range(p):
            left_node, right_node = self.split_on_j(splits,j)
            left_policy, left_score = left_node.find_best_policy_and_score(y)
            right_policy, right_score = right_node.find_best_policy_and_score(y)
            score = left_score + right_score
            if bestj is None or score > bestscore:
                bestj = j
                bestscore = score
        left_node, right_node = self.split_on_j(splits,bestj)
        return bestj, left_node, right_node

def get_splits(x):
    splits = []
    n = len(x)
    p = len(x[0])
    for j in range(p):
        left = []
        right = []
        for i in range(n):
            if x[i][j] == 0:
                left.append(i)
            else:
                right.append(i)
        splits.append((frozenset(left),frozenset(right)))
    return splits
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Learn an optimal policy tree from covariates and scores',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('data', help='file containing the CSV data (no commas in variables names allowed!)')
    parser.add_argument('--depth', default=2, type=int, help='(Maximum) depth of policy tree')
    parser.add_argument('--lazy', default=3, type=int, help='Laziness of unit policy determining constraints')
    parser.add_argument('--prefix', default=999999, type=int, help='Only use at most this many initial (unique) datapoints')
    parser.add_argument('--siblingcons', action="store_true", help='Whether to impose that neighbouring leaves have opposite policies')
    parser.add_argument('--usesos', action="store_true", help='Whether to use SOS1 constraints on branching nodes')
    parser.add_argument('--leftzerocons', action="store_true", help='Whether to impose that left leaves have policy 0 and right ones policy 1')
    parser.add_argument('--exclude', default="A,Y", help='Comma separated string of covariates to exclude')
    parser.add_argument('--scores', default="scores.DML,scores.DR,scores.CF", help='Comma separated string of fields to interpret as scores')
    parser.add_argument('--score', default="scores.DML", help='Score to use to build tree')
    
    args = parser.parse_args()

    data = np.genfromtxt(args.data,delimiter=',',skip_header=True)
    header = open(args.data).readline().rstrip()

    # use CSV reader to correctly parse header line
    
    for line in csv.reader([header],skipinitialspace=True):
        names = line
        break

    # construct x and y and get names for covariates
    
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

    # collapse identical covariate vectors
    
    x, indices  = np.unique(x,axis=0,return_inverse=True)
    newy = np.zeros(len(x))
    for i, idx in enumerate(indices):
        newy[idx] += y[i]
    y = newy
    
    # switch to lists
    
    x = x.tolist()[:args.prefix]
    y = y.tolist()[:args.prefix]

    splits = get_splits(x)
    n = len(x)
    p = len(x[0])

    root_indices = frozenset(range(n))
    
    root = Node(indices=root_indices,depth=0)
    root.tree_greedy(x,splits,y,args.depth)
    root.depth_first_enumerate(1,x_names)
    print(root.treestr())
    print('Score is {0}'.format(root.score()))

    root = Node(indices=root_indices,depth=0)
    root.tree_optimal(x,splits,y,args.depth)
    root.depth_first_enumerate(1,x_names)
    print(root.treestr())
    print('Score is {0}'.format(root.score()))

    model = Model()
    model.ModelSense = -1             # maximise objective (rather than the default of minimisation)
    model.Params.PreSOS1Encoding = 3
    root = Node(depth=0)
    root.maketree(args.depth)
    root.depth_first_enumerate(1,x_names)
    root.addmipvars(model,p,x_names)
    upvars = [model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i]) for i in range(n)]
    root.addmipcons(model,x,upvars,args.lazy,args.siblingcons,args.leftzerocons,args.usesos)
    model.optimize()
    root.apply_solution()
    print(root.treestr())
    print('Score is {0}'.format(model.objVal))
