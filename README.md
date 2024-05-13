This repo contains C code for constructing optimal policy trees from
covariate and reward data. The code here is used to create an R
package but one can create a standalone executable (if using Linux at
least).

To create an executable called "fpt" just do

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