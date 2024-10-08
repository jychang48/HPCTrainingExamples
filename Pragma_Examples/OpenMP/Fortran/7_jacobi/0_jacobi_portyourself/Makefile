EXEC = jacobi
default: ${EXEC}
all: ${EXEC}

ROCM_GPU ?= $(strip $(shell rocminfo |grep -m 1 -E gfx[^0]{1} | sed -e 's/ *Name: *//'))

FC1=$(notdir $(FC))

#ifeq ($(findstring amdflang,$(FC1)), amdflang)
#  OPENMP_FLAGS = -fopenmp
#  FREE_FORM_FLAG = -Mfreeform
ifeq ($(findstring -w amdflang-new,$(FC1)), amdflang-new)
  OPENMP_FLAGS = -fopenmp -fopenmp-version=52
  FREE_FORM_FLAG = -ffree-form
else ifeq ($(findstring -w flang,$(FC1)), flang)
  OPENMP_FLAGS = -fopenmp 
  FREE_FORM_FLAG = -Mfreeform
else ifeq ($(findstring flang-new,$(FC1)), flang-new)
  OPENMP_FLAGS = -fopenmp
  FREE_FORM_FLAG = -ffree-form
else ifeq ($(findstring gfortran,$(FC1)), gfortran)
  OPENMP_FLAGS = -fopenmp
  FREE_FORM_FLAG = -ffree-form
else ifeq ($(findstring ftn,$(FC1)), ftn)
  OPENMP_FLAGS = -fopenmp
endif

ifeq ($(findstring ftn,$(FC1)), ftn)
  FFLAGS = -G2 -O3 ${FREE_FORM_FLAG} ${OPENMP_FLAGS}
else ifeq ($(findstring amdflang-new,$(FC1)), amdflang-new)
  FFLAGS = -g -O3 -ffast-math ${FREE_FORM_FLAG} ${OPENMP_FLAGS}
else ifeq ($(findstring flang-new,$(FC1)), flang-new)
  FFLAGS = -g -O3 -ffast-math ${FREE_FORM_FLAG} ${OPENMP_FLAGS}
else
  FFLAGS = -g -O3 ${FREE_FORM_FLAG} -fstrict-aliasing ${OPENMP_FLAGS}
endif
ifeq (${FC1},gfortran-13)
   LDFLAGS = ${OPENMP_FLAGS} -fno-lto
else
   LDFLAGS = ${OPENMP_FLAGS}
endif

OBJS= \
	prec.o \
	mesh.o \
	input.o \
	norm.o \
	laplacian.o \
	boundary.o \
	update.o \
	jacobi.o \
	main.o \

boundary.o: prec.o mesh.o
input.o: prec.o mesh.o
jacobi.o: prec.o input.o mesh.o norm.o laplacian.o boundary.o update.o
laplacian.o: prec.o mesh.o
norm.o: prec.o mesh.o
update.o: prec.o mesh.o
main.o: mesh.o jacobi.o

%.o : %.f90 Makefile
	$(FC) $(FFLAGS) -c $< -o $@

${EXEC}: ${OBJS}
	$(FC) $(LDFLAGS) $^ -o $@

# Cleanup
clean:
	rm -f *.o ${EXEC} *.mod
