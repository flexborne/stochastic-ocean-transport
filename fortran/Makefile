FC ?= gfortran
FFLAGS ?= -O2 -Wall

all: transport_solver

transport_solver: config_module.o transport_solver.o
	$(FC) $(FFLAGS) -o $@ $^

config_module.o: config_module.f90
	$(FC) $(FFLAGS) -c $<

transport_solver.o: transport_solver.f90 config_module.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod transport_solver
