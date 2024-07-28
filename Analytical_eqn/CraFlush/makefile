# compiler
FC = gfortran

# compile flags
FCFLAGS = -I ./obj -J ./obj -g -fbounds-check -fbacktrace -W -Wunused -ffpe-trap=invalid,zero,overflow -ffpe-summary=none
# Using optimization (flag "-O") may give quite different results...
# link flags
FLFLAGS = -I ./obj -static

# "-I dir" Affects interpretation of the INCLUDE directive. This path is also used to
#   search for .mod files when previously compiled modules are required by a USE statement
# "-J dir" This option specifies where to put .mod files for compiled modules. It is
#   also added to the list of directories to searched by an USE statement.
#
# Note: hyphens not allowed in the "exe" name

go: go2 CRAFLUSH

go2:  # write a heading to the user
	@echo Usage is:    nmake  [CRAFLUSH]  [clean]

clean:
	del /q .\obj\*.*

CRAFLUSH: CraFlush.exe

CraFlush.exe : ./obj/CraFlush.o
	$(FC) $(FLFLAGS) -o CraFlush.exe -I ./obj ./obj/CraFlush.o

./obj/CraFlush.o : ./src/CraFlush.f90
	$(FC) $(FCFLAGS) -o ./obj/CraFlush.o -c ./src/CraFlush.f90
