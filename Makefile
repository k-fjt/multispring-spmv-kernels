FC = nvfortran
FFLAGCPU = -O3 -mp=multicore -Minfo=all -llapack
FFLAGGPU = -O3 -acc=gpu -gpu=cc90,cuda12.3,fastmath,loadcache:L1,lineinfo,deepcopy -Minfo=all -mp=multicore

# options
FLAG1 = -DNVEC=1 -DCRS10
FLAG2 = -DNVEC=1 -DCRS10
FLAG3 = -DNVEC=1

PROGRAM1 = kernelCRS_CPU_OPENMP.exe
PROGRAM2 = kernelCRS_GPU_OPENACC.exe
PROGRAM3 = kernelEBE_GPU_OPENACC.exe

SRCS = \
	main.F90 \
	kernelCRS.F \
	kernelEBE.F 

OBJS = $(SRCS:.F=.o)

.SUFFIXES: .o .F .F90

all: $(PROGRAM1) $(PROGRAM2) $(PROGRAM3) 

$(PROGRAM1):
	make clean
	$(FC) $(FLAG1) $(FFLAGCPU) -c $(SRCS)
	$(FC) $(FLAG1) $(FFLAGCPU) $(OBJS) -o $@

$(PROGRAM2):
	make clean
	$(FC) $(FLAG2) $(FFLAGGPU) -c $(SRCS)
	$(FC) $(FLAG2) $(FFLAGGPU) $(OBJS) -o $@

$(PROGRAM3):
	make clean
	$(FC) $(FLAG3) $(FFLAGGPU) -c $(SRCS)
	$(FC) $(FLAG3) $(FFLAGGPU) $(OBJS) -o $@

###################################################################################################
clean:
	rm -f *.o *.out *.lst *.s
clobber:
	rm -f *.o *.out *.lst *.s *.exe

