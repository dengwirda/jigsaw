# `JIGSAW: An unstrutured mesh generator`

<p align="center">
  <img src = "../master/img/bunny-TRIA3-1.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-2.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-3.png"> &nbsp
  <img src = "../master/img/bunny-TRIA4-3.png">
</p>

<a href="https://sites.google.com/site/dengwirda/jigsaw">`JIGSAW`</a> is a Delaunay-based unstructured mesh generator for two- and three-dimensional geometries. It is designed to generate high-quality triangular and tetrahedral meshes for planar, surface and volumetric problems. `JIGSAW` is based on a recently developed "restricted" Frontal-Delaunay algorithm -- a hybrid technique combining many of the best features of advancing-front and Delaunay-refinement type approaches.

`JIGSAW` is a stand-alone mesh generator written in C++, and is currently available for 64-bit `Windows` and `Linux` platforms. This package provides a basic command-line interface. A <a href="http://www.mathworks.com">`MATLAB`</a> // <a href="http://www.gnu.org/software/octave">`OCTAVE`</a> based scripting interface, including file I/O, mesh visualisation and post-processing facilities can be found <a href="https://github.com/dengwirda/jigsaw-matlab">here</a>.

# `Installation`

`JIGSAW` itself is a fully self-contained executable, without dependencies on third-party libraries or run-time packages. 

# `Starting Out`

`JIGSAW` is a command-line executable, and operates by reading and writing a set of input/output files. Run-time control for `JIGSAW` is enforced by instantiating the executable with a simple text-based configuration file (a so-called `JIG` file). Additional information on `JIGSAW`'s file-formats can be found <a href="http://sites.google.com/site/dengwirda/jigsaw/jigsaw-documentation">here</a>.

To run `JISAW`, first download and unzip the current <a href="https://github.com/dengwirda/jigsaw/archive/master.zip">repository</a>, then navigate to the installation directory and execute the following command-line entries:
````
On WIN-64 platforms:
\bin\WIN-64\jigsaw64r.exe example.jig

On GLX-64 platforms:
/bin/GLX-64/jigsaw64r     example.jig
````
Additional information, documentation, online tutorials and references are available <a href="http://sites.google.com/site/dengwirda/jigsaw">here</a>.


