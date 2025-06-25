# Fast Policy Tree

This repo provides R/C code for constructing optimal policy trees from
covariate and reward data. It aims to do the same job as the
`policytree` package but to do so more quickly.

# R package

Version 1.0 of the software (corresponding to commit `9ba0be2` in this repo) is [available on CRAN](https://cran.r-project.org/package=fastpolicytree).

## Using the R package

There is only one function: `fastpolicytree` which you can call
similarly to the `policy_tree` function from the `policytree`
package. As an example, the R script below finds optimal policy trees from
synthetic data with 500 datapoints, 30 binary covariates and rewards
for 2 actions. Depth 4 trees are found using policytree,
fastpolicytree and also
the [sparsepolicytree](https://github.com/beniaminogreen/sparsepolicytree) package.

``` r
library(fastpolicytree)
library(sparsepolicytree)
library(policytree)
Sys.setenv(RAYON_NUM_THREADS=1)

compare <- function(s, n, p, actions, depth, nvals)
{
        set.seed(s)
        X <- matrix(sample(nvals, n*p, replace = T), nrow=n, ncol=p)
        W <- sample(seq(1:actions), n, replace = TRUE)-1
        Y <- X[,1] + X[,2] * (W == 1) + X[,3] * (W == max(actions)) + runif(n, min=0)

        cf <- grf::causal_forest(X, Y, W)
        gamma <- double_robust_scores(cf)

        times <- list()
        for (method in 1:3)
        {
            if (method == 1 )
                time <- system.time(tree <- policy_tree(X, gamma, depth))
            else if ( method == 2 )
                time <- system.time(tree <- fastpolicytree(X, gamma, depth))
            else
                time <- system.time(tree <- sparse_policy_tree(X, gamma, depth))
            time <- as.vector(time)[3]
            times  <- c(times,time)
        }

        cat(s, n, p, actions, depth, times[[1]], times[[2]], times[[3]], "\n")
    }

nseeds  <- 5
n <- 500
p  <- 30
nactions  <- 2
depth  <- 4
nvals  <- 2

for (s in 1:nseeds)
    compare(s, n, p, nactions, depth, nvals)
```

Running this code on a 15 Gb RAM, 2.7GHz CPU laptop using Version 1.0
of fastpolicytree produces the following output:

```
1 500 30 2 4 137.829 0.861 278.163 
2 500 30 2 4 126.009 0.691 273.862 
3 500 30 2 4 126.13 0.788 268.526 
4 500 30 2 4 124.3 0.801 268.441 
5 500 30 2 4 124.057 0.823 268.486 
```

where the 6th column is the solving time for `policytree`, the 7th the
solving time for `fastpolicytree` and the 8th and final column the
solving time for `sparsepolicytree`. The trees will be the same for
each of the 3 methods. (This can be checked using the script
[test_getstats.R](R/test_getstats.R).)

The test above is unfair to `sparsepolicytree` since line
`Sys.setenv(RAYON_NUM_THREADS=1)` limits computation to a single CPU
and `sparsepolicytree` can exploit multiple CPUs. Changing this line
to `Sys.setenv(RAYON_NUM_THREADS=4)` produces the following output:

```
1 500 30 2 4 131.145 0.778 99.964 
2 500 30 2 4 134.879 0.696 99.654 
3 500 30 2 4 132.334 0.808 99.543 
4 500 30 2 4 133.471 0.825 99.507 
5 500 30 2 4 132.01 0.871 100.097 
```


## Installing the latest version

If you want to use the latest version (rather than version 1.0) and
you have a C compiler on your machine you can do:

``` r
devtools::install_github("jcussens/tailoring/fastpolicytree")
```

# Creating a standalone executable

It is possible to create a standalone executable (if using Linux at
least).

To create an executable called ``fpt`` just do

```
make
```

To get brief help on how to use just do

```
./fpt
```

To create a version with debugging turned on, do

```
make OPT=dbg
```

To remove compiled code do:

```
make clean
```

If you have doxygen installed then you can do

```
make doc
```

to create html and latex documentation.

# Funding

The development of `fastpolicytree` was supported by UK MRC project
  [Tailoring health policies to improve outcomes using machine
  learning, causal inference and operations research
  methods](https://gtr.ukri.org/projects?ref=MR%2FT04487X%2F1)

The following people worked on the "Tailoring ..." project

- [Noémi Kreif (PI)](https://sop.washington.edu/people/noemi-kreif/)
- [James Cussens](https://jcussens.github.io/)
- [Julia Hatamyar](https://www.york.ac.uk/che/staff/research/julia-hatamyar./)
- [Vishalie Shah](https://www.york.ac.uk/che/staff/students/vishalie-shah/)

# Licence, etc

All code in this repo:

- was written by [James Cussens](https://jcussens.github.io/). Please
  contact him at
  [james.cussens@bristol.ac.uk](mailto:james.cussens@bristol.ac.uk)
  with bug reports, questions, etc.
- is Copyright University of Bristol 2024
- uses GNU General Public Licence v3.0 (see LICENCE.txt)