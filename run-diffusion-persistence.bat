@set BIN=ph-compute\ScaleVariantTopo\Release\

@set DATA=GN
@set INPUT=demo\distmat\%DATA%-net

@set DIM=1
@set TAUMAX=0

@set NBEGIN=0
@set NEND=0
@set NORM=0

@set DST=demo\ph_%DATA%_norm_%NORM%
%BIN%DiffusionRunner_D64.exe --norm %NORM% --nthreads 16 -i %INPUT% -d %DIM% -t %TAUMAX% -o %DST% -b %NBEGIN% -e %NEND%
