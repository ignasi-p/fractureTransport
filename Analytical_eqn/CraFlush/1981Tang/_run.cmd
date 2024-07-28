@echo off
D:
cd D:\OneDrive\Documents\Geokemi\Modelling\PhreeqC\1D_transport\Analytical_eqn\craflush\1981Tang
cd

copy "craf_1E-11_alpha=0.in" "craf.in" >nul
..\craflush
copy "conc.out" "conc_1E-11_alpha=0.out" >nul
del craf.in > nul
del conc.out > nul

goto xit

copy craf_1E-10.in craf.in >nul
..\craflush
copy conc.out conc_1E-10.out >nul
del craf.in > nul
del conc.out > nul

copy craf_1E-14.in craf.in >nul
..\craflush
copy conc.out conc_1E-14.out >nul
del craf.in > nul
del conc.out > nul

:xit
pause
