## OpenPoroBEM

OpenPoroBEM is a Boundary Element Method (BEM) code for multiphase porous media subject to transient and dynamic loadings. The code can be used for dynamic (seismic site effects) and quasi-static analyses of bounded, unbounded, and half-space media. 

gfortran compiler:
```
$ sudo apt install gfortran-9
```

Compile from source code:
```
$ gfortran OpenPoroBEM.f90 -o OpenPoroBEM 
```

Run application:
```
$ ./ OpenPoroBEM 
```
 
