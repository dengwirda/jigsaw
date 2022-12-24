## `JIGSAW: An unstructured mesh generator`

<p align="middle">
  <img src = "../master/img/bunny-TRIA3-1.png" width="20%" hspace="0.25%">
  <img src = "../master/img/bunny-TRIA3-2.png" width="20%" hspace="0.25%">
  <img src = "../master/img/bunny-TRIA3-3.png" width="20%" hspace="0.25%">
  <img src = "../master/img/bunny-TRIA4-3.png" width="20%" hspace="0.25%">
</p>

`JIGSAW` is an unstructured mesh generator and tessellation library; designed to generate high-quality triangulations and polyhedral decompositions of general planar, surface and volumetric domains. `JIGSAW` includes refinement-based algorithms for the construction of new meshes, optimisation-driven techniques for the improvement of existing grids, as well as routines to assemble (restricted) Delaunay tessellations, Voronoi complexes and Power diagrams. 

This package provides the underlying `c++` source for `JIGSAW`; defining a basic command-line interface and a `c`-format `API`. Higher-level scripting interfaces, supporting a range of additional facilities for file I/O, mesh visualisation and post-processing operations are also available, including for <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="http://www.gnu.org/software/octave">`OCTAVE`</a> <a href="https://github.com/dengwirda/jigsaw-matlab">here</a> and for <a href="https://www.python.org/">`PYTHON`</a> <a href="https://github.com/dengwirda/jigsaw-python">here</a>.

`JIGSAW` is compiled and tested on various `64-bit` `Linux`, `Windows` and `MacOS` platforms using the `g++`, `clang++` and `msvc` compilers.

### `Code Structure`

`JIGSAW` is written as a `header-only` library in `c++`. Both a basic command-line interface and a `c`-format `API` are defined:

      JIGSAW::
      ├── src -- JIGSAW src code
      ├── inc -- JIGSAW header files (for libjigsaw)
      ├── bin -- JIGSAW's exe binaries live here
      ├── lib -- JIGSAW's lib binaries live here
      ├── geo -- geometry definitions and input data
      ├── out -- default folder for JIGSAW output
      └── uni -- unit tests and libjigsaw example programs

### `Getting Started`

The first step is to compile and configure the code! `JIGSAW` can either be built directly from src, or installed using the <a href="https://anaconda.org/conda-forge/jigsaw">`conda`</a> package manager.

### `Building from src`

The full `JIGSAW` src can be found in <a href="../master/src/">`../jigsaw/src/`</a>. It has been built using various `c++17` conforming versions of the `g++`, `clang++` and `msvc` compilers.

`JIGSAW` is a `header-only` package - the single main `jigsaw.cpp` file simply `#include`'s the rest of the library directly. `JIGSAW` does not currently dependent on any external packages or libraries.

`JIGSAW` consists of several pieces: `(a)` a set of command-line utilities that read and write mesh data from/to file, and `(b)` a shared library, accessible via a `c`-format `API`.

### `Using cmake`

`JIGSAW` can be built using the <a href="https://cmake.org/">`cmake`</a> utility. To build, follow the steps below:

    * Clone or download this repository.
    * Navigate to the root `../jigsaw/` directory.
    * Make a new temporary directory BUILD.
    * cd build
    * cmake .. -DCMAKE_BUILD_TYPE=BUILD_MODE
    * cmake --build . --config BUILD_MODE --target install EXTRAS
    * Delete the temporary BUILD directory.

This process will build a series of executables and shared libraries: `jigsaw` itself - the main command-line meshing utility, `tripod` - `JIGSAW`'s tessellation infrastructure, `marche` - a fast-marching solver designed to optimise mesh-spacing configurations, as well as `libjigsaw` - `JIGSAW`'s shared `API`. 

`BUILD_MODE` can be used to select different compiler configurations and should generally either be `Release` or `Debug`. `EXTRAS` can be used to pass additional compile-time arguments, for example `-- -j 4` will build in parallel on supported architectures.

See `example.jig` for documentation on calling the command-line executables, and the headers in <a href="../master/inc/">`../jigsaw/inc/`</a> for details on the `API`.

### `Using conda`

`JIGSAW` is also available as a `conda` environment. To install and use, follow the steps below:

    * Ensure you have conda installed. If not, consider miniconda as a lightweight option.
    * Add conda-forge as a channel: conda config --add channels conda-forge
    * Create a jigsaw environment: conda create -n jigsaw jigsaw

Each time you want to use `JIGSAW` simply activate the environment using: `conda activate jigsaw`

Once activated, the various `JIGSAW` command-line utilities will be available in your run path, `JIGSAW`'s shared library (`libjigsaw`) will be available in your library path and its include files in your include path.

### `CMD-line Examples`

After compiling the code, try running the following command-line example to get started:
````
On WIN platforms:

\bin\jigsaw.exe example.jig

On LNX platforms:

/bin/jigsaw     example.jig
````
In this example, a high-quality tetrahedral mesh is generated for the 'stanford-bunny' geometry and the result written to file. The input geometry is specified as a triangulated surface, and is read from `../jigsaw/geo/bunny.msh`. The volume and surface mesh outputs are written to `../jigsaw/out/bunny.msh`. See the `example.jig` text-file for a description of `JIGSAW`'s configuration options. 

A repository of additional surface models generated using `JIGSAW` can be found <a href="https://github.com/dengwirda/jigsaw-models">here</a>. A description of the `*.jig` and `*.msh` input file formats can be found in the <a href="https://github.com/dengwirda/jigsaw/wiki">wiki</a>.

### `libJIGSAW Scripts`

A set of unit-tests and `libjigsaw` example programs are contained in <a href="../master/uni/">`../jigsaw/uni/`</a>. The `JIGSAW-API` is documented via the header files in <a href="../master/inc/">`../jigsaw/inc/`</a>. 

The unit-tests can be built using the <a href="https://cmake.org/">`cmake`</a> utility. To build, follow the steps below:

    * Navigate to the `../jigsaw/uni/` directory.
    * Make a new temporary directory BUILD.
    * cd build
    * cmake .. -DCMAKE_BUILD_TYPE=BUILD_MODE
    * cmake --build . --config BUILD_MODE --target install EXTRAS
    * Delete the temporary BUILD directory.

This process will build the unit-tests as a series of executables in <a href="../master/uni/">`../jigsaw/uni/`</a>. `BUILD_MODE` is a compiler configuration flag: either `Release` or `Debug`. `EXTRAS` can be used to pass additional compile-time arguments.

### `Contributors`

1. [@dengwirda](github.com/dengwirda) is `JIGSAW`'s developer and maintainer --- this work was originally the focus of my PhD at the University of Sydney.
2. [@xylar](github.com/xylar) contributed the original `cmake` build system.
3. [@tunnellm](https://github.com/tunnellm) extended the original sequential optimisation algorithms to support thread-parallelism.

### `License`

This program may be freely redistributed under the condition that the copyright notices (including this entire header) are not removed, and no compensation is received through use of the software.  Private, research, and institutional use is free.  You may distribute modified versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS`. Distribution of this code as part of a commercial system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR`. (If you are not directly supplying this code to a customer, and you are instead telling them how they can obtain it for free, then you are not required to make any arrangement with me.) 

`DISCLAIMER`: Neither I nor `THE CONTRIBUTORS` warrant this code in any way whatsoever.  This code is provided "as-is" to be used at your own risk.

`THE CONTRIBUTORS` include:
(a) The University of Sydney
(b) The Massachusetts Institute of Technology
(c) Columbia University
(d) The National Aeronautics & Space Administration
(e) Los Alamos National Laboratory

### `References`

There are a number of publications that describe the algorithms used in `JIGSAW` in detail. If you make use of `JIGSAW` in your work, please consider including a reference to the following:

`[1]` - Darren Engwirda: Generalised primal-dual grids for unstructured co-volume schemes, J. Comp. Phys., 375, pp. 155-176, https://doi.org/10.1016/j.jcp.2018.07.025, 2018.

`[2]` - Darren Engwirda, Conforming Restricted Delaunay Mesh Generation for Piecewise Smooth Complexes, Procedia Engineering, 163, pp. 84-96, https://doi.org/10.1016/j.proeng.2016.11.024, 2016.

`[3]` - Darren Engwirda, Voronoi-based Point-placement for Three-dimensional Delaunay-refinement, Procedia Engineering, 124, pp. 330-342, http://dx.doi.org/10.1016/j.proeng.2015.10.143, 2015.

`[4]` - Darren Engwirda, David Ivers, Off-centre Steiner points for Delaunay-refinement on curved surfaces, Computer-Aided Design, 72, pp. 157-171, http://dx.doi.org/10.1016/j.cad.2015.10.007, 2016.

`[5]` - Darren Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, http://hdl.handle.net/2123/13148, 2014.


