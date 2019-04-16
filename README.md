## `JIGSAW: An unstructured mesh generator`

<p align="center">
  <img src = "../master/img/bunny-TRIA3-1.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-2.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-3.png"> &nbsp
  <img src = "../master/img/bunny-TRIA4-3.png">
</p>

`JIGSAW` is a computational library for unstructured mesh generation and tessellation; designed to generate high-quality triangulations and polyhedral decompositions of general planar, surface and volumetric domains. `JIGSAW` includes refinement-based algorithms for the construction of new meshes, optimisation-driven techniques for the improvement of existing grids, as well as routines to assemble (restricted) Delaunay tessellations, Voronoi complexes and Power diagrams.

This package provides the underlying `C++` source for `JIGSAW`; defining a basic command-line interface and a `C`-format `API`. A <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="http://www.gnu.org/software/octave">`OCTAVE`</a> based scripting interface, including a range of additional facilities for file I/O, mesh visualisation and post-processing operations can be found <a href="https://github.com/dengwirda/jigsaw-matlab">here</a>.

`JIGSAW` has been compiled and tested on various `64-bit` `Linux`, `Windows` and `MacOS` based platforms. 

## `Code Structure`

`JIGSAW` is written as a `header-only` library in `C++`. Both a basic command-line interface and a `C`-format `API` are defined:

      JIGSAW::
      ├── src -- JIGSAW src code
      ├── inc -- JIGSAW header files (for libjigsaw)
      ├── bin -- put JIGSAW exe binaries here
      ├── lib -- put JIGSAW lib binaries here
      ├── geo -- geometry definitions and input data
      ├── out -- default folder for JIGSAW output
      └── uni -- unit tests and libjigsaw example programs

## `Getting Started`

The first step is to compile and configure the code! `JIGSAW` can either be built directly from src, or installed using the <a href="https://anaconda.org/conda-forge/jigsaw">`conda`</a> package manager.

### `Building JIGSAW from src`

The full `JIGSAW` src can be found in <a href="../master/src/">`../jigsaw/src/`</a>.

`JIGSAW` is a `header-only` package - the single main `jigsaw.cpp` file simply `#include`'s the rest of the library directly. `JIGSAW` does not currently dependent on any external packages or libraries.

`JIGSAW` consists of several pieces: `(a)` a set of command-line utilities that read and write mesh data from/to file, and `(b)` a shared library, accessible via a `C`-format `API`.

#### `Using cmake`

`JIGSAW` can be built using the <a href="https://cmake.org/">`cmake`</a> utility. To build, follow the steps below:

    * Ensure you have the cmake utility installed.
    * Clone or download this repository.
    * Navigate to the root `../jigsaw/` directory.
    * Create a new temporary directory BUILD (to store the cmake build files).
    * Navigate into the temporary directory.
    * Execute: cmake -D CMAKE_BUILD_TYPE=BUILD_MODE ..
    * Execute: make
    * Execute: make install
    * Delete the temporary directory.

This process will build a series of executables and the shared library: `jigsaw` itself - the main command-line meshing utility, `tripod` - `JIGSAW`'s tessellation infrastructure, as well as `libjigsaw` - `JIGSAW`'s shared `API`. `BUILD_MODE` can be used to select different compiler configurations and should be either `RELEASE` or `DEBUG`. 

See `example.jig` for documentation on calling the command-line executables, and the headers in <a href="../master/inc/">`../jigsaw/inc/`</a> for details on the `API`.

#### `Via g++/llvm`

`JIGSAW` has been successfully built using various versions of the `g++` and `llvm` compilers. The build process is a simple one-liner (from <a href="../master/src/">`../jigsaw/src/`</a>):
````
g++ -std=c++11 -pedantic -Wall -O3 -flto -D NDEBUG
-D __cmd_jigsaw jigsaw.cpp -o ../bin/jigsaw
````
will build the main `jigsaw` command-line executable,
````
g++ -std=c++11 -pedantic -Wall -O3 -flto -D NDEBUG
-D __cmd_tripod jigsaw.cpp -o ../bin/tripod
````
will build the `tripod` command-line utility (`JIGSAW`'s tessellation infrastructure) and,
````
g++ -std=c++11 -pedantic -Wall -O3 -flto -fPIC -D NDEBUG
-D __lib_jigsaw jigsaw.cpp -shared -o ../lib/libjigsaw.so
````
will build `JIGSAW` as a shared library (`libjigsaw`).

### `Install JIGSAW via conda`

`JIGSAW` is also available as a `conda` environment. To install and use, follow the steps below:

	* Ensure you have conda installed. If not, consider miniconda as a lightweight option.
	* Add conda-forge as a channel: conda config --add channels conda-forge
	* Create a jigsaw environment: conda create -n jigsaw jigsaw

Each time you want to use `JIGSAW` simply activate the environment using: `conda activate jigsaw`

Once activated, the various `JIGSAW` command-line utilities will be available in your run path, `JIGSAW`'s shared library (`libjigsaw`) will be available in your library path and its include files in your include path.

## `A Simple Example`

After compiling the code, try running the following command-line example to get started:
````
On WIN platforms:

\bin\jigsaw.exe example.jig

On LNX platforms:

/bin/jigsaw     example.jig
````
In this example, a high-quality tetrahedral mesh is generated for the 'stanford-bunny' geometry and the result written to file. The input geometry is specified as a triangulated surface, and is read from `../jigsaw/geo/bunny.msh`. The volume and surface mesh outputs are written to `../jigsaw/out/bunny.msh`. See the `example.jig` text-file for a description of `JIGSAW`'s configuration options. 

A repository of additional surface models generated using `JIGSAW` can be found <a href="https://github.com/dengwirda/jigsaw-models">here</a>.

## `API Example Code`

A set of unit-tests and `libjigsaw` example programs are contained in <a href="../master/uni/">`../jigsaw/uni/`</a>. The `JIGSAW-API` is documented via the header files in <a href="../master/inc/">`../jigsaw/inc/`</a>. 

The unit-tests can be built using the <a href="https://cmake.org/">`cmake`</a> utility. To build, follow the steps below:

    * Navigate to the `../jigsaw/uni/` directory.
    * Create a new temporary directory BUILD (to store the cmake build files).
    * Navigate into the temporary directory.
    * Execute: cmake -D CMAKE_BUILD_TYPE=BUILD_MODE ..
    * Execute: make
    * Execute: make install
    * Delete the temporary directory.

This process will build the unit-tests as a series of executables in <a href="../master/uni/">`../jigsaw/uni/`</a>. `BUILD_MODE` is a compiler configuration flag: either `RELEASE` or `DEBUG`.

## `License`

This program may be freely redistributed under the condition that the copyright notices (including this entire header) are not removed, and no compensation is received through use of the software.  Private, research, and institutional use is free.  You may distribute modified versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS`. Distribution of this code as part of a commercial system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR`. (If you are not directly supplying this code to a customer, and you are instead telling them how they can obtain it for free, then you are not required to make any arrangement with me.) 

`DISCLAIMER`:  Neither I nor: Columbia University, the Massachusetts Institute of Technology, the University of Sydney, nor the National Aeronautics and Space Administration warrant this code in any way whatsoever.  This code is provided "as-is" to be used at your own risk.

## `References`

There are a number of publications that describe the algorithms used in `JIGSAW` in detail. If you make use of `JIGSAW` in your work, please consider including a reference to the following:

`[1]` - Darren Engwirda: Generalised primal-dual grids for unstructured co-volume schemes, J. Comp. Phys., 375, pp. 155-176, https://doi.org/10.1016/j.jcp.2018.07.025, 2018.

`[2]` - Darren Engwirda, Conforming Restricted Delaunay Mesh Generation for Piecewise Smooth Complexes, Procedia Engineering, 163, pp. 84-96, https://doi.org/10.1016/j.proeng.2016.11.024, 2016.

`[3]` - Darren Engwirda, Voronoi-based Point-placement for Three-dimensional Delaunay-refinement, Procedia Engineering, 124, pp. 330-342, http://dx.doi.org/10.1016/j.proeng.2015.10.143, 2015.

`[4]` - Darren Engwirda, David Ivers, Off-centre Steiner points for Delaunay-refinement on curved surfaces, Computer-Aided Design, 72, pp. 157-171, http://dx.doi.org/10.1016/j.cad.2015.10.007, 2016.

`[5]` - Darren Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, http://hdl.handle.net/2123/13148, 2014.


