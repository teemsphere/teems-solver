all: all-redirect

BUILD_DIR ?= /opt/teems-solver
PETSC_DIR = $(BUILD_DIR)/lib/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC        = mpicc
MPICH_DIR = $(BUILD_DIR)/lib/mpi
FC        = $(MPICH_DIR)/bin/mpif90 -fopenmp
F77       = $(MPICH_DIR)/bin/mpif77 -fopenmp
OPT      ?= -Ofast
F77FLAGS  = -I. -I/usr/local/include
CFLAGS    = -I. -I${PETSC_DIR}/native/include -I ${PETSC_DIR}/include -I${PETSC_DIR}/src/mat/impls/aij/seq -fopenmp
FFLAGS    = -fallow-argument-mismatch
LDFLAGS   = -fopenmp

LIBS      = ${PETSC_KSP_LIB} $(BUILD_DIR)/lib/ma48/src/.libs/libma48.so $(BUILD_DIR)/lib/ma51/src/.libs/libma51.so

all-redirect: teems-solver
OBJECTS = main.o globals.o cmf_io.o str_util.o tab_parse.o formula.o jacobian.o block_order.o block_solve.o solve_drivers.o solve_rk.o probe.o hsl_kernels.o hsl_mp48d.o ddeps.o hsl_mc66d.o ddeps90.o hsl_mc79i.o hsl_mc79i_ciface.o
teems-solver: $(OBJECTS)
	${CLINKER} -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(OPT) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) $(OPT) -c $<
%.o: %.f
	$(F77) $(F77FLAGS) $(OPT) -c $<

.PHONY: clean
# PETSc's clean:: (clean-legacy) already removes *.o and *.mod. That
# also deletes hsl_mp01.mod, which is a staged input from the mp48
# library build rather than a product of this makefile, so restore it.
clean::
	-rm -f *.o
	-rm -f teems-solver
	@test -f hsl_mp01.mod || cp $(BUILD_DIR)/lib/mp48/src/hsl_mp01.mod . 2>/dev/null || true

# dependencies
hsl_mp48d.o: ddeps.o
hsl_mc66d.o: ddeps90.o
hsl_mc79i_ciface.o: hsl_mc79i.o
hsl_kernels.o: hsl_mc66d.o hsl_mp48d.o
