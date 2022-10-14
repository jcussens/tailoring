# tailoring

C and Python code for building policy trees. All code written by James Cussens

For the C code to compile do the following (from Linux command line):
    `gcc -DNDEBUG -O2 opttree.c main.c`
which will produce an executable file `a.out`
    
To run `a.out` do:
    ` ./a.out IFLS.txt 2 2`

The penultimate argument is the number of actions.
The final argument is the desired tree depth.
    
Here IFLS.txt is a file whose first line is a header giving the names of the covariates, followed by the names of the actions. Each subsequent line has the covariate values for a unit, followed by reward values for that unit for each action. Values are separated by spaces.

Here is an example of producing a depth 2 tree with 10,622 units, 64 binary covariates and 2 actions. (Note that it takes under a second on James's laptop.)

(base) uw20605@IT079795:~/repos/tailoring$ time ./a.out IFLS.txt 2 2
Actions: 0: "scores.DML" 1: "scores.DR" 
node = 0x55913d0cb480
covariate = "past_ind_mis"
value = 1
reward = 69.3649
action_id = 0
left_child = 0x55913d0cb4c0
right_child = 0x55913d0cb580

node = 0x55913d0cb4c0
covariate = "province.f18"
value = 1
reward = 32.0115
action_id = 0
left_child = 0x55913d0cb500
right_child = 0x55913d0cb540

node = 0x55913d0cb500
value = 0
reward = 0
action_id = 1
left_child = (nil)
right_child = (nil)

node = 0x55913d0cb540
value = 0
reward = 32.0115
action_id = 0
left_child = (nil)
right_child = (nil)

node = 0x55913d0cb580
covariate = "region.f3-Jawa"
value = 1
reward = 37.3533
action_id = 0
left_child = 0x55913d0cb5c0
right_child = 0x55913d0cb600

node = 0x55913d0cb5c0
value = 0
reward = 0
action_id = 1
left_child = (nil)
right_child = (nil)

node = 0x55913d0cb600
value = 0
reward = 37.3533
action_id = 0
left_child = (nil)
right_child = (nil)


real	0m0.890s
user	0m0.865s
sys	0m0.020s

This repo also contains the R script "pt.R" which can be used to run
policytree on input file with the same format as that accepted by the C code. Here is an example of finding a depth 2 tree on the same data as just above:

(base) uw20605@IT079795:~/repos/tailoring$ time Rscript pt.R IFLS.txt 2 2
Warning in policy_tree(x, gammas, depth) :
  The number of covariates exceeds 50. Consider reducing the dimensionality before running policy_tree, by for example using only the Xj's with the highest variable importance (`grf::variable_importance` - the runtime of exact tree search scales with ncol(X)^depth, see the documentation for details).
policy_tree object 
Tree depth:  2 
Actions:  1: scores.DML 2: scores.DR 
Variable splits: 
(1) split_variable: past_ind_mis  split_value: 0 
  (2) split_variable: province.f18  split_value: 0 
    (4) * action: 2 
    (5) * action: 1 
  (3) split_variable: region.f3.Jawa  split_value: 0 
    (6) * action: 2 
    (7) * action: 1 

real	5m26.860s
user	5m26.449s
sys	0m0.210s

We get the same tree as with the C version, but it takes about 327 seconds rather than less than 1 second.
