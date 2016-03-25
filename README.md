# `JIGSAW: An unstrutured mesh generator`

<p align="center">
  <img src = "../master/img/bunny-TRIA3-1.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-2.png"> &nbsp
  <img src = "../master/img/bunny-TRIA3-3.png"> &nbsp
  <img src = "../master/img/bunny-TRIA4-3.png">
</p>

<a href="https://sites.google.com/site/dengwirda/jigsaw">`JIGSAW`</a> is a Delaunay-based unstructured mesh generator for two- and three-dimensional geometries. It is designed to generate high-quality triangular and tetrahedral meshes for planar, surface and volumetric problems. `JIGSAW` is based on a recently developed "restricted" Frontal-Delaunay algorithm -- a hybrid technique combining many of the best features of advancing-front and Delaunay-refinement type approaches.

`JIGSAW` is a stand-alone mesh generator written in C++, currently available for 64-bit `Windows` and `Linux` platforms. This package provides a basic command-line interface. A <a href="http://www.mathworks.com">`MATLAB`</a> // <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> based scripting interface, including file I/O, mesh visualisation and post-processing facilities can be found <a href="">here<\a>.

# `Installation`

`JIGSAW` itself is a fully self-contained executable, without dependencies on third-party libraries or run-time packages. 

# `Starting Out`

After downloading and unzipping the current <a href="../master/archive/master.zip">repository</a>, navigate to the installation directory and run the following example problem:
````
meshdemo(1); % build surface-meshes
meshdemo(2); % build volume-meshes
meshdemo(3); % preserve "sharp-features" in piecewise smooth domains
meshdemo(4); % build planar-meshes -- impose topological constraints
meshdemo(5); % build planar meshes -- explore mesh-size controls
````
Additional information, documentation, online tutorials and references are available <a href="https://sites.google.com/site/dengwirda/jigsaw">here</a>.


