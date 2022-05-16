from gurobipy import *
import argparse
import numpy as np
import csv
import sys

class Node:
    '''Represents nodes of a policy tree.

    The root node of a tree implicitly represents the entire tree.

    If a non-leaf node then the covariate defining a split may be defined
    or may not. If a leaf node then the policy for that leaf may be defined
    or may not. So this class is used during the process of finding a good (perhaps optimal)
    policy tree and also for representing particular policy trees.
    '''
    def __init__(self,covariate=None,covariate_name=None,name=None,left=None,right=None,mipvars=None,
                 mother=None,indices=None,depth=None,policy=None):
        '''Initialise a policy tree node

        Args: 

         covariate (int/None): If None, this indicates that no
          'splitting' covariate has been chosen for this node. 
          Otherwise the splitting covariate specified either by index.
         covariate_name (str): If None, this indicates that no
          'splitting' covariate has been chosen for this node. 
          Otherwise the splitting covariate specified by name.
         name (str): The name of the node
         left (Node/None): For a non-leaf, the left daughter of the node, else None.
         right (Node/None): For a non-leaf, the right daughter of the node, else None.
         mipvars (list): The set of MIP vars representing possible choices of covariate (if any).
         mother (Node/None): For a non-root node, the mother of the node, else None.
         indices (frozenset): The set of (indices of) datapoints that reach this node if this can be determined, else None.
         depth (int): The depth of the node. Root nodes have depth 0.
         policy (int): The policy if the node is a leaf node with an decided policy, else None.
        '''
        self._covariate = covariate
        self._covariate_name = covariate_name
        self._name = name
        self._left = left
        self._right = right
        self._mipvars = mipvars
        self._mother = mother
        self._indices = indices
        self._depth = depth
        self._policy = policy

    def is_leaf_node(self):
        '''Return whether the node is a leaf node

        Returns:
         bool: Whether the node is a leaf node
        '''
        assert self._left is None or self._right is not None
        assert self._right is None or self._left is not None
        return self._left is None and self._right is None


    def is_branching_node(self):
        '''Return whether the node is a branching node

        Returns:
         bool: Whether the node is a branching node
        '''
        return not self.is_leaf_node()

    
    def mother(self):
        '''Return the mother of the node (which is None for root node)

        Returns:
         Node/None: Mother of the node
        '''
        return self._mother

    def children(self):
        '''Return the children of the node if there are any, or None

        Returns:
         tuple/None: (left child, right child) / None
        '''
        if self._left is None:
            return None
        else:
            return self._left, self._right

    def sibling(self):
        '''Return the sibling of the node, or None if it has no sibling

        Returns:
         Node: The sibling of the node, or None if it has no sibling
        '''
        if self._mother is None:
            return None
        left, right = self._mother.children()
        if left is self:
            return right
        else:
            return left

    
    def left_child(self):
        '''Return the left child of a node if it has one, else None

        Returns:
         Node: Left child of a node
        '''
        return self._left

    def leaf_policy_mipvar(self):
        '''Return the MIP variable specifying the policy for the leaf node

        Raises:
         ValueError: If the node is not a leaf node

        Returns:
         (gurobipy.Var): The MIP variable specifying the policy for the leaf node
        '''
        if self._left is not None:
            raise ValueError("Not a leaf node, so has no policy")
        return self._mipvars[0]
    
    def is_left_child(self):
        '''Returns whether the node is the left child of some other node

        Returns:
         bool: Whether the node is the left child of some other node
        '''
        return self._mother is not None and self._mother.left_child() is self

    def set_covariate(self,covariate):
        '''Set the (splitting) covariate for a node

        Args:
         covariate (int): The covariate specified by index
        '''
        self._covariate = covariate

    def set_children(self,left_child,right_child):
        '''Set the children for a node

        Args:
         left_child (Node): The left child
         right_child (Node): The right child
        '''
        self._left = left_child
        self._right = right_child

    def set_policy(self,policy):
        '''Make a node a leaf node with given policy

        Args:
         policy (int): The policy (0 or 1)
        '''
        self._covariate = None
        self._left = None
        self._right = None
        self._policy = policy

    def score(self):
        '''Return the score for the tree rooted at the node

        Returns:
         float: The score for the tree rooted at the node
        '''
        if self._policy == 0:
            return 0
        elif self._policy == 1:
            return sum([y[i] for i in self._indices])
        else:
            return self._left.score() + self._right.score()

    def treestr(self):
        '''Returns a string representing the (sub-)tree for which the node is the root.

        Returns:
         str: String representing the (sub-)tree for which the node is the root.    
        '''
        res = str(self)
        if self._left is not None:
            res += self._left.treestr()
        if self._right is not None:
            res += self._right.treestr()
        return res

    def __str__(self):
        '''Returns a string representing the node

        Returns:
         str: String representing the node
        '''
        res = '{0}\n'.format(self._name)
        if self.is_leaf_node():
            res += 'policy   : {0}\n'. format(self._policy)
        else:
            res += 'covariate: {1}({0})\n'. format(self._covariate,self._covariate_name)
            res += 'left     : {0}\n'. format(self._left._name)
            res += 'right    : {0}\n'. format(self._right._name)
        return res + '\n'

    def get_path_from_root(self):
        '''Returns the path from the root node to the node

        Path is a list of (node,decision) pairs. For example, if a node
        was reached by two left branches followed by a right branch the path
        would be of the form [(root,0),(node1,0),(node2,1)]

        Returns:
         list: Path from the root node to the node
        '''
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
        '''
        Apply a solution to the MIP model to fix a particular policy tree

        Raises:
         AttributeError: If no MIP solution available.
        '''
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
            if self.leaf_policy_mipvar().X > 0.5:
                self._policy = 1
            else:
                self._policy = 0
            
    def addmipvars(self,model,x_names):
        '''Add MIP variables to the MIP model

        Args:
         model (gurobipy.Model): The MIP model
         x_names (iter): The names of the covariates
        '''
        if self.is_branching_node():
            self._mipvars = []
            for x_name in x_names:
                v = model.addVar(name="a#{0}#{1}".format(self._name,x_name),vtype=GRB.BINARY)
                v.BranchPriority = 10-self._depth
                self._mipvars.append(v)
            nvars_added_left = self._left.addmipvars(model,x_names)
            nvars_added_right = self._right.addmipvars(model,x_names)
            nvars_added = len(x_names) + nvars_added_left + nvars_added_right
        else:
            v = model.addVar(name="np#{0}".format(self._name),vtype=GRB.BINARY)
            v.BranchPriority = 100
            self._mipvars = [v]
            nvars_added = 1
        return nvars_added
    
    def addmipcons(self,model,x,upvars,lazy,siblingcons,leftzerocons,usesos):
        '''Add constraints to the MIP model

        Args:
         model (gurobipy.Model): The MIP model
         x (np.ndarray): Matrix of covariate values, one row for each unit, one column
          for each covariate
         upvars (list): upvars[i] is the MIP variable representing the policy for unit i
         lazy (int): The 'lazy' attribute for constraints connecting covariate vectors to 
          leaves and policies. (There are many such constraints so they have to be introduced
          lazily as solving progresses. This value controls how readily the constraint is pulled into the model.)
         siblingcons (bool): Whether to add the constraint that sibling leaves always have opposite
          policies.
         leftzerocons (bool): Whether to fix left leaf polices to 0 and right leaf policies to 1.
         usesos (bool): Whether to add SOS-1 constraints for each branching node. 
        '''
        if self.is_branching_node():
            # exactly one splitting attribute for each branching node
            model.addConstr(quicksum(self._mipvars) == 1)
            if usesos:
                model.addSOS(GRB.SOS_TYPE1, self._mipvars)
            self._left.addmipcons(model,x,upvars,lazy,siblingcons,leftzerocons,usesos)
            self._right.addmipcons(model,x,upvars,lazy,siblingcons,leftzerocons,usesos)
        else:
            leaf_policy_mipvar = self.leaf_policy_mipvar()
            if self._policy is not None:
                model.addConstr(leaf_policy_mipvar == self._policy)

            if self.is_left_child():
                if leftzerocons:
                    # just fix policy to 0 and set sibling's policy to 1
                    model.addConstr(leaf_policy_mipvar == 0)
                    model.addConstr(self.sibling().leaf_policy_mipvar() == 1)
                elif siblingcons:
                    # if a left child state that it has opposite policy from sibling
                    model.addConstr(leaf_policy_mipvar + self.sibling().leaf_policy_mipvar() == 1)
                
            nodesdecisions = self.get_path_from_root()
            # for each unit that reaches this leaf ensure policy for unit is the same as policy for leaf
            for i, xi in enumerate(x):
                lexpr = LinExpr()
                for node, decision in nodesdecisions:
                    avars = [v for j, v in enumerate(node._mipvars) if xi[j] == decision]
                    lexpr.addConstant(1)
                    lexpr.addTerms([-1]*len(avars),avars)
                # lexpr is 0 iff xi reaches this leaf, otherwise it is at least 1,
                # and the following constraints are redundant
                unit_policy_mipvar = upvars[i]
                if self._policy is None:
                    cons = model.addConstr(lexpr + leaf_policy_mipvar - unit_policy_mipvar >= 0)
                    cons.Lazy = lazy
                    cons = model.addConstr(lexpr - leaf_policy_mipvar + unit_policy_mipvar >= 0)
                    cons.Lazy = lazy
                elif self._policy == 0:
                    cons = model.addConstr(lexpr + leaf_policy_mipvar - unit_policy_mipvar >= 0)
                    cons.Lazy = lazy
                else:
                    cons = model.addConstr(lexpr - leaf_policy_mipvar + unit_policy_mipvar >= 0)
                    cons.Lazy = lazy
                    

    def maketree(self,depth):
        '''Grow a tree of the given depth below the given node
        
        Creates tree structure, but does not fix splitting attributes, leaf policies.
        Does not add any MIP variables or constraints.
        If depth is 0 no tree is created.

        Args:
         depth (int): Depth of tree to grow
        '''
        if depth == 0:
            return
        left_node = Node(mother=self,depth=self._depth+1)
        right_node = Node(mother=self,depth=self._depth+1)
        self.set_children(left_node,right_node)
        left_node.maketree(depth-1)
        right_node.maketree(depth-1)
            
    def depth_first_enumerate(self,idx,x_names):
        '''Give nodes integer labels in a depth-first manner.

        Also assigns covariate names

        Args:
         idx (int): label for the root node
         x_names (list) : x_names[j] is the name of the jth covariate 

        Returns:
         int: the next available integer label
        '''
        self._name = str(idx)
        if self._covariate is not None:
            self._covariate_name = x_names[self._covariate]
        next = idx + 1
        if self._left is not None:
            next = self._left.depth_first_enumerate(next,x_names)
        if self._right is not None:
            next = self._right.depth_first_enumerate(next,x_names)
        return next

    def split_on_j(self,splits,j):
        '''Split (indices of) the data on a covariate and return resulting daughter nodes.

        Args:
         splits (iter): jth element is (go_left_indices,go_right_indices) where go_left_indices is the set of indices of datapoints
          that go left (have value 0 for covariate j), and go_right_indices similarly for having value 1
         j (int): Covariate to split on

        Returns:
         (tuple): A left daughter node and right daughter node which corresponding sets of data indices associated with each.
        '''
        left_indices =  splits[j][0] & self._indices
        right_indices = splits[j][1] & self._indices
        return Node(mother=self,indices=left_indices,depth=self._depth+1), Node(mother=self,indices=right_indices,depth=self._depth+1)
        
    def find_best_policy_and_score(self,y):
        '''Find best policy and score of that policy if node were to become a leaf node

        Args:
         y (iter): Vector of response values
        '''
        score_for_1 = sum([y[i] for i in self._indices])
        # we are minimising
        if score_for_1 < 0:
            return 1, score_for_1
        else:
            return 0, 0

    def tree_greedy(self,x,splits,y,treedepth):
        '''Find a policy tree for the given subset of the data and depth
        using the greedy algorithm
        
        Args:
         x (np.ndarray): Matrix of covariate values, one row for each unit in the data subset, one column
          for each covariate
         splits (iter): jth element is (go_left_indices,go_right_indices) where go_left_indices is the set of indices of datapoints
          that go left (have value 0 for covariate j), and go_right_indices similarly for having value 1
         y (iter): Vector of response values
         treedepth (int): Desired tree depth
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
        '''Find an optimal policy tree for the given subset of the data and depth

        The tree is not returned. It is specified by setting attributes of nodes.

        Args:
         x (np.ndarray): Matrix of covariate values, one row for each unit in the data subset, one column
          for each covariate
         splits (iter): jth element is (go_left_indices,go_right_indices) where go_left_indices is the set of indices of datapoints
          that go left (have value 0 for covariate j), and go_right_indices similarly for having value 1
         y (iter): Vector of response values
         treedepth (int): Desired tree depth

        Returns:
         float: Score of optimal tree
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
                # we are minimising
                if bestj is None or score < best_score:
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

        Args:
         x (np.ndarray): Matrix of covariate values, one row for each unit in the data subset, one column
          for each covariate
         splits (iter): jth element is (go_left_indices,go_right_indices) where go_left_indices is the set of indices of datapoints
          that go left (have value 0 for covariate j), and go_right_indices similarly for having value 1
         y (iter): Vector of response values

        Returns:
         tuple: (best_covariate, left_node, right_node)
        '''
        bestj = None
        p = len(x[0])  # assume non-empty data
        for j in range(p):
            left_node, right_node = self.split_on_j(splits,j)
            left_policy, left_score = left_node.find_best_policy_and_score(y)
            right_policy, right_score = right_node.find_best_policy_and_score(y)
            score = left_score + right_score
            if bestj is None or score < bestscore:
                bestj = j
                bestscore = score
        left_node, right_node = self.split_on_j(splits,bestj)
        return bestj, left_node, right_node

                
def get_splits_j(x):
    '''For each covariate get indices of datapoints
    that go left and those that go right

    Args:
     x (np.ndarray): Matrix of covariate values, one row for each unit, one column
          for each covariate

    Returns:
     list: jth element is a pair (tuple) of sets, the set of ('left') data indices that have value
      0 for covariate j, followed by the set of ('right') data indices that have value
      1 for covariate j.
    '''
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
    parser.add_argument('--greedyno', action="store_true", help='Whether *not* to run the greedy algorithm')
    parser.add_argument('--nonmipno', action="store_true", help='Whether *not* to run the exact non-MIP algorithm')
    parser.add_argument('--mipno', action="store_true", help='Whether *not* to run the exact MIP algorithm')
    parser.add_argument('--verbose', action="store_true", help='Whether to be verbose')
    
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

    n = len(x)
    p = len(x[0])

    if args.verbose:
        print('{0} units and {1} covariates.'.format(n,p))
    
    splits_j = get_splits_j(x)
    root_indices = frozenset(range(n))

    if not args.greedyno:
        root = Node(indices=root_indices,depth=0)
        root.tree_greedy(x,splits_j,y,args.depth)
        root.depth_first_enumerate(1,x_names)
        print(root.treestr())
        print('Score is {0}'.format(root.score()))

    if not args.nonmipno:
        root = Node(indices=root_indices,depth=0)
        root.tree_optimal(x,splits_j,y,args.depth)
        root.depth_first_enumerate(1,x_names)
        print(root.treestr())
        print('Score is {0}'.format(root.score()))

    if not args.mipno:
        model = Model()
        model.Params.PreSOS1Encoding = 3
        root = Node(depth=0)
        root.maketree(args.depth)
        root.depth_first_enumerate(1,x_names)
        nvars_added = root.addmipvars(model,x_names)
        upvars = [model.addVar(name="up#{0}".format(i),vtype=GRB.BINARY,obj=y[i]) for i in range(n)]
        nvars_added += n
        if args.verbose:
            print('{0} MIP vars in model, of which {1} are for unit policy'.format(nvars_added,n))
        root.addmipcons(model,x,upvars,args.lazy,args.siblingcons,args.leftzerocons,args.usesos)
        model.optimize()
        root.apply_solution()
        print(root.treestr())
        print('Score is {0}'.format(model.objVal))


# POSSIBLE FUTURE CODE TO HELP CONSTRUCT BETTER MIP FORMULATION

    #splits_i = get_splits_i(x)


    # splitsi = get_splits_i(x)
    # #print(splitsi)
    # relations(splitsi)

    #best_split_exact(frozenset(range(n)),p,splits_i)

    #best_split(frozenset(range(n)),p,splits_j,splits_i)

    #find_best_j(x,p,frozenset(range(n)),splits_j)


# def best_split_exact(indices,p,splits_i):
#     model = Model()
#     # VARIABLES
#     b = []
#     for j in range(p):
#         b.append(model.addVar(name='b#{0}'.format(j),vtype=GRB.BINARY))
#     l = {}
#     r = {}
#     u = {}
#     for i in indices:
#         l[i] = model.addVar(name='l#{0}'.format(i),vtype=GRB.BINARY)
#         r[i] = model.addVar(name='r#{0}'.format(i),vtype=GRB.BINARY)
#         u[i] = model.addVar(name='u#{0}'.format(i),vtype=GRB.BINARY,obj=1)
#     # CONSTRAINTS
#     size = int(p/2)
#     for i in indices:
#         model.addConstr(l[i] + r[i] + u[i] == 1)
#     model.addConstr(quicksum(b) == size)
#     for i in indices:
#         li = splits_i[i][0]
#         ri = splits_i[i][1]
#         maxl = min(size,len(li))
#         maxr = min(size,len(ri))
#         model.addConstr(quicksum([b[j] for j in li]) <= maxl - maxl*r[i])
#         model.addConstr(quicksum([b[j] for j in ri]) <= maxr - maxr*l[i])
#         #model.addConstr((r[i]==1) >> (quicksum([b[j] for j in li]) == 0))
#         #model.addConstr((l[i]==1) >> (quicksum([b[j] for j in ri]) == 0))
#     model.optimize()
    
        
# def get_splits_i(x):
#     '''For each datapoint x[i] return those covariates for which x[i][j]=0
#     and those for which x[i][j]=1
#     '''
#     splits = []
#     n = len(x)
#     p = len(x[0])
#     for i in range(n):
#         xi = x[i]
#         left = []
#         right = []
#         for j in range(p):
#             if xi[j] == 0:
#                 left.append(j)
#             else:
#                 right.append(j)
#         splits.append((frozenset(left),frozenset(right)))
#     return splits

# def get_initial(splits_j):
#     maxsplitsize = 0
#     for j, (left,right) in enumerate(splits_j):
#         splitsize = max(len(left),len(right))
#         if splitsize > maxsplitsize:
#             maxsplitsize = splitsize
#             best_j = j
#     return best_j

# def best_split(indices,p,splits_j,splits_i):
#     remaining = set(range(p))
#     existing = set()
#     best_candidate = get_initial(splits_j)
#     existing.add(best_candidate)
#     remaining.remove(best_candidate)
#     while len(existing) <= p/2:
#         best_candidate, best_determined = greedy_choose_j(indices,splits_i,existing,remaining)
#         existing.add(best_candidate)
#         remaining.remove(best_candidate)
#         print(existing,best_determined)
        
# def get_determined(indices,new_set,splits_i):
#     determined = 0
#     for i in indices:
#         if new_set <= splits_i[i][0] or new_set <= splits_i[i][1]:
#             determined += 1
#     return determined
            
# def greedy_choose_j(indices,splits_i,existing,remaining):
#     best_determined = 0
#     for candidate in remaining:
#         determined = get_determined(indices,existing | frozenset([candidate]),splits_i)
#         if determined > best_determined:
#             best_determined  = determined
#             best_candidate = candidate
#     return best_candidate, best_determined

# def find_best_j(x,p,indices,splits_j):
#     if len(indices) < 2:
#         return
#     goal_size = len(indices)/2
#     best_size = goal_size
#     best_j = None
#     for j in range(p):
#         left = indices & splits_j[j][0]
#         this_size = abs(len(left) - goal_size)
#         if this_size < best_size:
#             best_size = this_size
#             best_j = j
#             if best_size < 1:
#                 break
#     best_left = indices & splits_j[best_j][0]
#     best_right = indices & splits_j[best_j][1]
#     print('Splitting {0} datapoints on {1} into a left size of {2} and a right size of {3}'.format(len(indices),best_j,len(best_left),len(best_right)))
#     find_best_j(x,p,best_left,splits_j)
#     find_best_j(x,p,best_right,splits_j)

            
# def relations(splits):
#     for i1, (left1,right1) in enumerate(splits):
#         for tmp, (left2,right2) in enumerate(splits[i1+1:]):
#             i2 = tmp + i1 + 1
#             if right1 <= right2:
#                 print(i1,i2,'{0} <= {1}'.format(i1,i2))
#             elif right2 <= right1:
#                 print(i1,i2,'{0} >= {1}'.format(i1,i2))
#             elif len(right1 & right2) == 0:
#                 print(i1,i2,'{0} & {1} = 0'.format(i1,i2))
#             else:
#                 for tmp2, (left3,right3) in enumerate(splits[i2+1:]):
#                     if (right1 <= right3 or right3 <= right1 or right2 <= right3 or right3 <= right2 or
#                         len(right1&right3)==0 or len(right2&right3)==0):
#                         continue
#                     i3 = tmp2 + i2 + 1
#                     if right1 & right2 <= right3:
#                         print(i1,i2,i3,'{0} & {1} <= {2}'.format(i1,i2,i3))
#                     if right1 & right3 <= right2:
#                         print(i1,i2,i3,'{0} & {1} <= {2}'.format(i1,i3,i2))
#                     if right2 & right3 <= right1:
#                         print(i1,i2,i3,'{0} & {1} <= {2}'.format(i2,i3,i1))
#                     if len(right1 & right2 & right3) == 0:
#                         print(i1,i2,i3,'{0} & {1} & {2} = 0'.format(i1,i2,i3))
