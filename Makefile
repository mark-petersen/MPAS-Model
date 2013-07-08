.SUFFIXES: .F .o

OBJS = mpas_rbf_interpolation.o \
       mpas_vector_reconstruction.o \
       mpas_spline_interpolation.o \
       mpas_tracer_advection_helpers.o \
       mpas_tracer_advection_mono.o \
       mpas_tracer_advection_std.o \
	   mpas_geometry_utils.o

all: operators

operators: $(OBJS)
	ar -ru libops.a $(OBJS)

mpas_vector_reconstruction.o: mpas_rbf_interpolation.o
mpas_rbf_interpolation.o:
mpas_spline_interpolation:
mpas_tracer_advection_helpers.o: mpas_geometry_utils.o
mpas_tracer_advection_mono.o: mpas_tracer_advection_helpers.o
mpas_tracer_advection_std.o: mpas_tracer_advection_helpers.o 
mpas_geometry_utils.o:

clean:
	$(RM) *.o *.mod *.f90 libops.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
