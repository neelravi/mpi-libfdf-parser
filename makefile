#
#  Makefile for FDF library and FDF example.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: check_flags


#
# Architecture
#
FDF_ARCH=XLF 64bits PARALLEL

#
# Compiling, linking flags
#
FC= mpiifort
FFLAGS= -O0 -fpp -g -debug -traceback # -qstrict -qtune=ppc970 -qarch=ppc970 -q32
LDFLAGS=-O0 #-qstrict -qtune=ppc970 -qarch=ppc970 -q32


#
# Library utils
#
#AR=
RANLIB= echo

#
# MPI flags
#
MPIFLAGS= $(MPI_PATH) $(MPI_INCLUDE) $(MPI_LIBS)

MPI_PATH=
MPI_INCLUDE=
MPI_LIBS=

#
# Macro definitions
#
# DEBUG:    Full debugging information in library
# _MPI_:    Runs FDF under a MPI environment
# CLUSTER:  FDF on a Cluster (non-shared filesystem)
# BLOCKING: Blocking reading on a huge MPI execution
#
DEFS= $(DEFS_DEBUG) $(DEFS_MPI) $(DEFS_IO)

DEFS_DEBUG= -DDEBUG
DEFS_MPI= -D_MPI_
DEFS_IO= -DCLUSTER # or #-WF,-DBLOCKING


COMPILE = $(FC) $(FFLAGS) $(DEFS) -c
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location
# of the arch.make file.
#
#
#ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Add items to the INCFLAGS as set by the arch.make file
#
#INCFLAGS:=-I$(VPATH) $(INCFLAGS)    # For VPATH operation
#
# This is needed to pick up top-level module files, such as mpi_siesta
#
#LOCAL_INCFLAGS=-I../
#INCFLAGS:=$(INCFLAGS) $(LOCAL_INCFLAGS)
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libfdf.a


SOURCES=iso_fortran_env.F90 io_fdf.F90 parse.F90 utils.F90 prec.F90 fdf.F90

#SOURCES=iso_fortran_env.F90 io_fdf.F90 parse.F90 utils.F90 prec.F90 mpi_write.F90 fdf.F90


%.o %.mod %.smod: %.F90
	$(COMPILE) -o $*.o $<
	@touch $@

main: $(subst .f90,.o,$(SOURCES))
	$(FC) -o $@ $+


#$(subst .F90,.o,$(SOURCES))

default: module
module: libfdf.a
	cp libfdf.a ..
	@cp -p *.*d ..
#
lib: iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o
	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a \
                iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o
	-$(RANLIB) libfdf.a

#lib:  mpi_write.o iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o
#	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a \
#                mpi_write.o iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o
#	-$(RANLIB) libfdf.a



#check_flags:
#	@echo "In fdf, INCFLAGS is: $(INCFLAGS)"
#
fdf.o: io_fdf.o parse.o prec.o utils.o
io_fdf.o: utils.o prec.o
parse.o: utils.o prec.o
utils.o: prec.o

#
# Tests for FDF library (Sequential, MPI, IO)
#
sample: sample.F90 libfdf.a
	$(FC) $(FFLAGS) $(LDFLAGS) $(MPIFLAGS) $(DEFS) -o $@ $^

#
# Cleaning
#
clean:
	rm -f *.o *.mod libfdf.a
	rm -f sample sample_mpi sample_cluster sample_blocking
	rm -f fdf.debug* sample.out* sample.fdf.*
	rm -f *.i90
