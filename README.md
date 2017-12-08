## `JIGSAW: An unstrutured mesh generator`

<p align="center">
  <img src = "../master/img/bunny-TRIA3-1.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-2.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-3.png"> &nbsp
  <img src = "../master/img/bunny-TRIA4-3.png">
</p>

`JIGSAW` is a computational library for unstructured mesh generation; designed to generate high-quality triangulations and polyhedral decompositions of general planar, surface and volumetric domains. `JIGSAW` includes both `refinement`-based algorithms for the construction of new meshes, as well as `optimisation`-driven techniques for the improvement of existing grids.

This package provides the underlying `C++` source for `JIGSAW`; defining a basic command-line interface and a `C`-format `API`. A <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="http://www.gnu.org/software/octave">`OCTAVE`</a> based scripting interface, including a range of additiona facilities for file I/O, mesh visualisation and post-processing operations can be found <a href="https://github.com/dengwirda/jigsaw-matlab">here</a>.

`JIGSAW` has been compiled and tested on various `64-bit` `Linux` , `Windows` and `Mac` based platforms. 


## `Code Structure`

`JIGSAW` is written in `C++` as a `header-only` library. Support is provided for a basic command-line interface and a `C`-format `API`.

      JIGASW::
      ├── src -- JIGSAW src code
      ├── inc -- JIGSAW header files (for libjigsaw)
      ├── bin -- put JIGSAW exe binaries here
      ├── lib -- put JIGSAW lib binaries here
      ├── geo -- geometry definitions and input data
      ├── out -- default folder for JIGSAW output
      └── uni -- unit tests and libjigsaw example programs

## `Getting Started`

The first step is to compile the code! The `JIGSAW` src can be found in <a href="../master/jigsaw/src/">`../jigsaw/src/`</a>.

`JIGSAW` is a `header-only` package - there is only the single main `jigsaw.cpp` file that simply `#include`'s the rest of the library as headers. The resulting build process should be fairly straight-forward as a result. `JIGSAW` does not currently dependent on any external packages or libraries.

#### `On Linux/Mac`

`JIGSAW` has been successfully built using various versions of the `g++` and `llvm` compilers. Since the build process is a simple one-liner, there's no `make` script - instead:

	g++ -std=c++11 -pedantic -Wall -s -O3 -flto -D NDEBUG -I libcpp -static-libstdc++ 
	jigsaw.cpp -o jigsaw64r
	
can be used to build a `JIGSAW` executable, while:

	g++ -std=c++11 -pedantic -Wall -O3 -flto -fPIC -D NDEBUG -I libcpp -static-libstdc++ 
	jigsaw.cpp -shared -o libjigsaw64r.so

can be used to build a `JIGSAW` shared library. See the headers in <a href="../master/jigsaw/inc/">`../jigsaw/inc/`</a> for details on the `API`. The `#define __lib_jigsaw` directive in `jigsaw.cpp` toggles the source between executable and shared-library modes.

#### `On Windows`

`JIGSAW` has been successfully built using various versions of the `msvc` compiler. I do not provide a sample `msvc` project, but the following steps can be used to create one:

	* Create a new, empty MSVC project.
	* Import the jigsaw.cpp file, this contains the main() entry-point.
	* Modify the MSVC project settings to include the "../src/" and "../src/libcpp/" directories.

#### `Folder Structure`

Once you have built the `JIGSAW` binaries, place them in the appropriate sub-folders in`../jigsaw/bin/` and/or `../jigsaw/lib/` directories, so that they can be found by the unit tests in `../jigsaw/uni/`.

## `A Simple Example`

After compiling the code, run the following command-line example to confirm:
````
On WIN-64 platforms:

\bin\WIN-64\jigsaw64r.exe example.jig

On LNX-64 platforms:

/bin/LNX-64/jigsaw64r     example.jig
````
In this example, a high-quality tetrahedral mesh is generated for the 'stanford-bunny' geometry and the result is written to file. The input geometry is specified as a triangulated surface, and is read from `../geo/bunny.msh`. The volume and surface mesh outputs are written to `../out/bunny.msh`. Edit `example.jig` for a description of `JIGSAW`'s configuration options. 

Additional information, documentation, online tutorials and references are available <a href="http://sites.google.com/site/dengwirda/jigsaw">here</a>. A repository of 3D surface models generated using `JIGSAW` can be found <a href="https://github.com/dengwirda/jigsaw-models">here</a>.

## `Attribution!`

If you make use of `JIGSAW` please make reference to the following. The algorithmic developments behind `JIGSAW` have been the subject of a number of publications, originally stemming from my PhD research at the University of Sydney:

`[1]` - Darren Engwirda: Generalised primal-dual grids for unstructured co-volume schemes, under review, https://arxiv.org/abs/1712.02657, 2017.

`[2]` - Darren Engwirda, Conforming Restricted Delaunay Mesh Generation for Piecewise Smooth Complexes, Procedia Engineering, Volume 163, Pages 84-96, ISSN 1877-7058, https://doi.org/10.1016/j.proeng.2016.11.024, 2016.

`[3]` - Darren Engwirda, Voronoi-based Point-placement for Three-dimensional Delaunay-refinement, Procedia Engineering, Volume 124, Pages 330-342, ISSN 1877-7058, http://dx.doi.org/10.1016/j.proeng.2015.10.143, 2015.

`[4]` - Darren Engwirda, David Ivers, Off-centre Steiner points for Delaunay-refinement on curved surfaces, Computer-Aided Design, Volume 72, Pages 157-171, ISSN 0010-4485, http://dx.doi.org/10.1016/j.cad.2015.10.007, 2016.

`[5]` - Darren Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, http://hdl.handle.net/2123/13148, 2014.

