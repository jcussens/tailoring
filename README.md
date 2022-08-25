# tailoring

C and Python code for building policy trees. All code written by James Cussens

For the C code to compile do:
    `gcc -DNDEBUG opttree.c main.c`
which will produce an executable file `a.out`
    
To run `a.out` do:
    ` ./a.out IFLS.txt 2`
    
The final argument is the desired tree depth.
    
Here IFLS.txt is a file whose first 3 lines are 3 numbers, e.g.:

    10622
    64
    2
    
where line 1 is the number of units, line 2 is the number of covariates and line 3 is the number of actions. 
line 4 and onwards contain the data separated by white space. on each data line we have the covariate values followed by the reward values for each action
