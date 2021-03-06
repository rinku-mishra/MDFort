# MDFort (Molecular Dynamics in Fortran)

[![build](https://github.com/rinku-mishra/MDFort/actions/workflows/main.yml/badge.svg)](https://github.com/rinku-mishra/MDFort/actions/workflows/main.yml)

A OpenMp parallel Fortran code to simulate plasma particles using Molecular Dynamics Algorithm.

## Problem
<!--Rayleigh Problem = gas between 2 plates ([Alexander & Garcia, 1997](https://doi.org/10.1063/1.168619)) -->

## Contributors
- [Rinku Mishra](https://github.com/rinku-mishra), CPP-IPR, India. [@arra4723](https://twitter.com/arra4723)
- [Rupak Mukherjee](https://github.com/RupakMukherjee), PPPL, USA.

Installation
------------
#### Prerequisites
1. [GNU Make](https://www.gnu.org/software/make/)
2. [gfortran](https://gcc.gnu.org/fortran/)
3. [git](https://git-scm.com/)

#### Procedure
First make a clone of the master branch using the following command
```shell
git clone https://github.com/rinku-mishra/MDFort.git
```
Then enter inside the *MDFort* directory
```shell
cd MDFort
```
Now complile and built the *tara2d* code
```shell
make subsystems
``` 

```shell
make all
``` 
Usage
-----
Upon successful compilation, run the code using following command
```shell
make run
```
