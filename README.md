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

On LNX-64 platforms:

/bin/LNX-64/jigsaw64r     example.jig
````
In this example, a high-quality tetrahedral mesh is generated for the stanford-bunny geometry and the result is written to file. The input geometry is specified as a triangulated surface, and is read from `../geo/bunny.msh`. The volume and surface mesh outputs are written to `../out/bunny.msh`. Edit `example.jig` for a description of `JIGSAW`'s configuration options. 

Additional information, documentation, online tutorials and references are available <a href="http://sites.google.com/site/dengwirda/jigsaw">here</a>.

# `Attribution!`

If you make use of `JIGSAW` please reference appropriately. The algorithmic developments behind `JIGSAW` have been the subject of a number of publications, beginning with my PhD research at the University of Sydney:

`[1]` - Darren Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014, http://hdl.handle.net/2123/13148.

`[2]` - Darren Engwirda, David Ivers, Off-centre Steiner points for Delaunay-refinement on curved surfaces, Computer-Aided Design, Volume 72, March 2016, Pages 157-171, ISSN 0010-4485, http://dx.doi.org/10.1016/j.cad.2015.10.007.

`[3]` - Darren Engwirda, Voronoi-based Point-placement for Three-dimensional Delaunay-refinement, Procedia Engineering, Volume 124, 2015, Pages 330-342, ISSN 1877-7058, http://dx.doi.org/10.1016/j.proeng.2015.10.143. 


