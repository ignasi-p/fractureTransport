@echo off
rem echo Minimalist GNU for Windows with Fortran compiler
rem echo Minimalist GNU for Windows-64 with Fortran compiler  

echo Minimalist GNU for Windows with Fortran compiler
set PATH=C:\bin\mingw64\bin;%PATH%
rem echo %PATH%
echo working directory:
cd

echo.
mingw32-make CRAFLUSH

echo.
pause
