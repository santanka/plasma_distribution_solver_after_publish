FC      = gfortran
FCFLAGS = -O4
LDFLAGS = -fopenmp
TARGET  = a.out

OBJ  = constant_parameter.o \
       constant_in_the_simulation.o \
	   reference_results_setting.o \
       main_variables.o \
       subroutines.o \
       main.o \

all: $(TARGET)

.SUFFIXES : .f90

.f90.o:
	$(FC) -c $< $(FCFLAGS) $(LDFLAGS)

$(TARGET): $(OBJ)
	$(FC) $(OBJ) $(FCFLAGS) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.exe core