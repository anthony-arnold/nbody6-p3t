# NBODY6+P3T

This is a version of the [NBODY6](https://people.ast.cam.ac.uk/~sverre/web/pages/nbody.htm) program by Sverre Aarseth, modified to use the P3T algorithm described in [Iwasawa et al. (2015)](https://arxiv.org/abs/1506.04553)

The code forks from the NBODY6 source and is modified to use the [CMake](https://cmake.org/) build system.

## Papers

If you use this code, please cite the following papers:

 - Arnold, Baumgardt & Wang, [Accelerating NBODY6 with a graphics processing unit-enabled Particle-Particle Particle-Tree scheme](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2075A/abstract), MNRAS, 509, 2, 2022.
 - Arnold & Baumgardt, [Direct N-body simulations of NGC 6397 and its tidal tails](https://ui.adsabs.harvard.edu/abs/2025arXiv250111806A/abstract), MNRAS, 2025.

## Structure

The original structure of the code has been retained. Where modifications were required to a Fortran subroutine, a *.tree.f file was created to hold the modified version. Additionally, the [code](https://github.com/treecode/bonsai) for the GPU-enabled tree code [Bonsai](https://arxiv.org/abs/1204.2280) is included in the GPU2/lib directory.

## Building

To build everything, create a new build directory and run `cmake` from this directory with the GPU2 directory as the target and then build the targets. For example:

    $ git clone git@github.com:anthony-arnold/nbody6-p3t
    $ mkdir build
    $ cd build
    $ cmake ../nbody6-p3t/GPU2
    $ make -j4

CMake will generate two targets, `nbody7.gpu` and `nbody7.tree`. The former target refers to the original NBODY6 program and the latter is the new program.

## Input

The new program takes the same input as the original NBODY6 with a small change. After the initial input which is common to both programs, the following variables are read: `NNBOPT` and `RBUFF`. The first variable is unused and can be left as `0`. The second variable is the R_buff value required by the P3T algorithm.

Additionally, the opening angle value for the BH-tree algorithm is needed and is set in the environment variable `THETA`.
