@echo off
rem gfortran 6.0

SET fp=
SET fort=gfortran
SET opt= -std=f2008 -O2

set cl=%fort%%opt%
echo %cl%

%cl% ^
%fp%classCell.f90 ^
%fp%classSpatialInfo.f90 ^
%fp%classFluidModel.f90 ^
%fp%classTimeIterator.f90 ^
%fp%classOutputData.f90 ^
%fp%VOF2D.f90 ^
-static -o VOF2D

del *.mod *.obj