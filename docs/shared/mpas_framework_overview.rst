.. role:: raw-latex(raw)
   :format: latex
..

MPAS Framework Overview
=======================

The MPAS Framework provides the foundation for a generalized geophysical
fluid dynamics model on unstructured spherical and planar meshes. On top
of the framework, implementations specific to the modeling of a
particular physical system (e.g., land ice, ocean) are created as MPAS
*cores*. To date, MPAS cores for atmosphere, ocean, shallow
water, sea ice, and land ice have been implemented. The MPAS design
philosophy is to leverage the efforts of developers from the various
MPAS cores to provide common framework functionality with minimal
effort, allowing MPAS core developers to focus on development of the
physics and features relevant to their application.

The framework code includes shared modules for fundamental model
operation. Significant capabilities include:

-  *Description of model data types.* MPAS uses a handful of fundamental
   Fortran derived types for basic model functionality. Core-specific
   model variables are handled through custom groupings of model fields
   called *pools*, for which custom accessor routines exist.
   Core-specific variables are easily defined in XML syntax in a
   *Registry*, and the framework parses the Registry, defines variables,
   and allocates memory as needed.

-  *Description of the mesh specification.* MPAS requires 36 fields to
   fully describe the mesh used in a simulation. These include the
   position, area, orientation, and connectivity of all cells, edges,
   and vertices in the mesh. The mesh specification can flexibly
   describe both spherical and planar meshes. More details are provided
   in the next section.

-  *Distributed memory parallelization and domain decomposition.* The
   MPAS Framework provides needed routines for exchanging information
   between processors in a parallel environment using Message Passing
   Interface (MPI). This includes halo updates, global reductions, and
   global broadcasts. MPAS also supports decomposing multiple domain
   blocks on each processor to, for example, optimize model performance
   by minimizing transfer of data from disk to memory. Shared memory
   parallelization through OpenMP is also supported, but the
   implementation is left up to each core.

-  *Parallel input and output capabilities.* MPAS performs parallel
   input and output of data from and to disk through the commonly used
   libraries of NetCDF, Parallel NetCDF (pnetcdf), and Parallel
   Input/Output (PIO). The Registry
   definitions control which fields can be input and/or output, and a
   framework *streams* functionality provides easy run-time
   configuration of what fields are to be written to what file name and
   at what frequency through an XML streams file. The MPAS framework
   includes additional functionality specific to providing a flexible
   model restart capability.

-  *Advanced timekeeping.* MPAS uses a customized version of the
   timekeeping functionality of the Earth System Modeling Framework
   (ESMF), which includes a robust set of time and calendar tools used
   by many Earth System Models (ESMs). This allows explicit definition
   of model epochs in terms of years, months, days, hours, minutes,
   seconds, and fractional seconds and can be set to three different
   calendar types: Gregorian, Gregorian no leap, and 360 day. This
   flexibility helps enable multi-scale physics and simplifies coupling
   to ESMs. To manage the complex date/time types that ensue, MPAS
   framework provides routines for arithmetic of time intervals and the
   definition of alarm objects for handling events (e.g., when to write
   output, when the simulation should end).

-  *Run-time configurable control of model options.* Model options are
   configured through *namelist* files that use standard Fortran
   namelist file format, and input/output are configured through
   *streams* files that use XML format. Both are completely adjustable
   at run time.

-  *Online, run-time analysis framework.* A system for defining analysis
   of model states during run time, reducing the need for
   post-processing and model output.

Additionally, a number of shared operators exist to perform common
operations on model data. These include geometric operations (e.g.,
length, area, and angle operations on the sphere or the plane),
interpolation (linear, barycentric, Wachspress, radial basis functions,
spline), vector and tensor operations (e.g., cross products,
divergence), and vector reconstruction (e.g., interpolating from cell
edges to cell centers). Most operators work on both spherical and planar
meshes.
