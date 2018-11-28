.. role:: math(raw)
   :format: html latex
..

Grid Description
================

This chapter provides a brief introduction to the common types of grids
used in the MPAS framework.

The MPAS grid system requires the definition of seven elements. These
seven elements are composed of two types of *cells*, two types of
*lines*, and three types of *points*. These elements are depicted in
Figure [figure:variablePosition] and defined in Table
[table:variablePosition]. These elements can be defined on either the
plane or the surface of the sphere. The two types of cells form two
meshes, a primal mesh composed of Voronoi regions and a dual mesh
composed of Delaunay triangles. Each corner of a primal mesh cell is
uniquely associated with the “center” of a dual mesh cell and vice
versa. So we define the two mesh as either a primal mesh (composed of
cells :math:`P_i`) or a dual mesh (composed of cells :math:`D_v`). The
center of any primal mesh cell, :math:`P_i`, is denoted by
:math:`{\bf x}_i` and the center of any the dual mesh cell, :math:`D_v`,
is denoted by :math:`{\bf x}_v`. The boundary of a given primal mesh
cell :math:`P_i` is composed of the set of lines that connect the
:math:`{\bf x}_v` locations of associated dual mesh cells :math:`D_v`.
Similarly, the boundary of a given dual mesh cell :math:`D_v` is
composed of the set of lines that connect the :math:`{\bf x}_i`
locations of the associated primal mesh cells :math:`P_i`.

As shown in Figure [figure:variablePosition], a line segment that
connects two primal mesh cell centers is uniquely associated with a line
segment that connects two dual mesh cell centers. We assume that these
two line segments cross and the point of intersection is labeled as
:math:`{\bf x}_e`. In addition, we assume that these two line segments
are orthogonal as indicated in Figure [figure:variablePosition]. Each
:math:`{\bf x}_e` is associated with two distances: :math:`d_e` measures
the distance between the primal mesh cells sharing :math:`{\bf x}_e` and
:math:`l_e` measures the distance between the dual mesh cells sharing
:math:`{\bf x}_e`.

Since the two line segments crossing at :math:`{\bf x}_e` are
orthogonal, these line segments form a convenient local coordinate
system for each edge. At each :math:`{\bf x}_e` location a unit vector
:math:`{\bf n}_e` is defined to be parallel to the line connecting
primal mesh cells. A second unit vector :math:`{\bf t}_e` is defined
such that :math:`{\bf t}_e = {\bf k} \times {\bf n}_e`.

In addition to these seven element types, we require the definition of
*sets of elements*. In all, eight different types of sets are required
and these are defined and explained in Table [table:gridConnectivity]
and Figure [figure:gridConnectivity]. The notation is always of the form
of, for example, :math:`i \in CE(e)`, where the LHS indicates the type
of element to be gathered (cells) based on the RHS relation to another
type of element (edges).

Table [table:gridFileName] provides the names of all *elements* and all
*sets of elements* as used in the MPAS framework. Elements appear twice
in the table when described in the grid file in more than one way, e.g.
points are described with both cartesian and latitude/longitude
coordinates. An “ncdump -h” of any MPAS grid, output or restart file
will contain all variable names shown in second column of Table
[table:gridFileName].

+---------------------+----------------+------------------------------------------------------------+
| :math:`Element`     | :math:`Type`   | :math:`Definition`                                         |
+=====================+================+============================================================+
| :math:`{\bf x}_i`   | point          | location of center of primal-mesh cells                    |
+---------------------+----------------+------------------------------------------------------------+
| :math:`{\bf x}_v`   | point          | location of center of dual-mesh cells                      |
+---------------------+----------------+------------------------------------------------------------+
| :math:`{\bf x}_e`   | point          | location of edge points where velocity is defined          |
+---------------------+----------------+------------------------------------------------------------+
| :math:`d_{e}`       | line segment   | distance between neighboring :math:`{\bf x}_i` locations   |
+---------------------+----------------+------------------------------------------------------------+
| :math:`l_{e}`       | line segment   | distance between neighboring :math:`{\bf x}_v` locations   |
+---------------------+----------------+------------------------------------------------------------+
| :math:`P_i`         | cell           | a cell on the primal-mesh                                  |
+---------------------+----------------+------------------------------------------------------------+
| :math:`D_v`         | cell           | a cell on the dual-mesh                                    |
+---------------------+----------------+------------------------------------------------------------+

Table: Definition of elements used to build the MPAS grid.

+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`Syntax`           | :math:`ouptut`                                                                       |    |
+==========================+======================================================================================+====+
| :math:`e \in EC(i) `     | set of edges that define the boundary of :math:`P_i`.                                |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`e \in EV(v) `     | set of edges that define the boundary of :math:`D_v`.                                |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`i \in CE(e) `     | two primal-mesh cells that share edge :math:`e`.                                     |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`i \in CV(v) `     | set of primal-mesh cells that form the vertices of dual mesh cell :math:`D_v`.       |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`v\in VE(e) `      | the two dual-mesh cells that share edge :math:`e`.                                   |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`v \in VI(i) `     | the set of dual-mesh cells that form the vertices of primal-mesh cell :math:`P_i`.   |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`e \in ECP(e)`     | edges of cell pair meeting at edge :math:`e`.                                        |    |
+--------------------------+--------------------------------------------------------------------------------------+----+
| :math:`e \in EVC(v,i)`   | edge pair associated with vertex :math:`v` and mesh cell :math:`i`.                  |    |
+--------------------------+--------------------------------------------------------------------------------------+----+

Table: Definition of element groups used to reference connections in the
MPAS grid. Examples are provided in Figure [figure:gridConnectivity].

+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`Element`        | :math:`Name`      | :math:`Size`         | :math:`Comment`                                |
+========================+===================+======================+================================================+
| :math:`{\bf x}_i`      | {x,y,z}Cell       | nCells               | cartesian location of :math:`{\bf x}_i`        |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`{\bf x}_i`      | {lon,lat}Cell     | nCells               | longitude and latitude of :math:`{\bf x}_i`    |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`{\bf x}_v`      | {x,y,z}Vertex     | nVertices            | cartesian location of :math:`{\bf x}_v`        |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`{\bf x}_v`      | {lon,lat}Vertex   | nVertices            | longitude and latitude of :math:`{\bf x}_v`    |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`{\bf x}_e`      | {x,y,z}Edge       | nEdges               | cartesian location of :math:`{\bf x}_e`        |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`{\bf x}_e`      | {lon,lat}Edge     | nEdges               | longitude and latitude of :math:`{\bf x}_e`    |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`d_{e}`          | dcEdge            | nEdges               | distance between :math:`{\bf x}_i` locations   |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`l_{e}`          | dvEdge            | nEdges               | distance between :math:`{\bf x}_v` locations   |
+------------------------+-------------------+----------------------+------------------------------------------------+
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`e \in EC(i) `   | edgesOnCell       | (nEdgesMax,nCells)   | edges that define :math:`P_i`.                 |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`e \in EV(v) `   | edgesOnVertex     | (3,nCells)           | edges that define :math:`D_v`.                 |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`i \in CE(e) `   | cellsOnEdge       | (2,nEdges)           | primal-mesh cells that share edge :math:`e`.   |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`i \in CV(v) `   | cellsOnVertex     | (3,nVertices)        | primal-mesh cells that define :math:`D_v`.     |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`v\in VE(e) `    | verticesOnEdge    | (2,nEdges)           | dual-mesh cells that share edge :math:`e`.     |
+------------------------+-------------------+----------------------+------------------------------------------------+
| :math:`v \in VI(i) `   | verticesOnCell    | (nEdgesMax,nCells)   | vertices that define :math:`P_i`.              |
+------------------------+-------------------+----------------------+------------------------------------------------+

Table: Variable names used to describe a MPAS grid.

| |Definition of elements used to build the MPAS grid. Also see Table
  [table:variablePosition].|

| |Definition of element groups used to reference connections in the
  MPAS grid. Also see Table [table:gridConnectivity].|

.. |Definition of elements used to build the MPAS grid. Also see Table [table:variablePosition].| image:: ./shared/figures/variablePosition.pdf
   :width: 16.00000cm
.. |Definition of element groups used to reference connections in the MPAS grid. Also see Table [table:gridConnectivity].| image:: ./shared/figures/gridConnectivity.pdf
   :width: 16.00000cm
