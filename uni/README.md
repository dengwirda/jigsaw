## `UNIT-TESTS: various benchmarks for JIGSAW`

A set of unit tests designed to exercise various `JIGSAW`, `TRIPOD` and `MARCHE` functionality:

    * TEST2D_X.c: 2D Cartesian test cases, with X = {A, B, ..., etc}.
    * TEST3D_X.c: 3D Cartesian test cases, with X = {A, B, ..., etc}.
    * TEST2S_X.c: 2D Spherical test cases, with X = {A, B, ..., etc}.

    * TESTKZ_A.c: Mesh a simple domain.
    * TESTKZ_B.c: Mesh with user-defined h(x) input (a structured 'grid').
    * TESTKZ_C.c: Mesh with user-defined h(x) input (unstructured 'mesh').
    * TESTKZ_D.c: Mesh 'multi-part' geometry with tags.
    * TESTKZ_E.c: Mesh with imposed initial conditions.
    * TESTKZ_F.c: Build 'restricted' DTs, using TRIPOD.
    * TESTKZ_G.c: Build conventional DTs, using TRIPOD.
    * TESTKZ_H.c: Solve: |dh/dx| <= G, using MARCHE (a structured 'grid').
    * TESTKZ_I.c: Solve: |dh/dx| <= G, using MARCHE (unstructured 'mesh'). 

See `JIGSAW`'s base `README` for build information via <a href="https://cmake.org/">`cmake`</a>.


