@set SRC=%1
@set BIN=release\

@set DATA=GN
@set INPUT=..\distmat\%DATA%-net

@set DIM=1
@set TAUMAX=0

@set NBEGIN=0
@set NEND=0

@set NPROCS=4
@set NORM=0

@set DST=%SRC%ph_%DATA%_norm_%NORM%
%BIN%DiffusionRunner_D64.exe --norm %NORM% --nprocs %NPROCS% -i %INPUT% -d %DIM% -t %TAUMAX% -o %DST% --nbegin %NBEGIN% --nend %NEND%
