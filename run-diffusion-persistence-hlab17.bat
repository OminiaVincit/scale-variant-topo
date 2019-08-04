@set BIN=ph-compute\ScaleVariantTopo\Release\

@set DATA=REDDIT-BINARY
@set EXP_PATH=F:\tran\Research\ScaleVariant\exp_20190516\%DATA%
@set LABEL=norm_1_avg_1

@set INPUT=%EXP_PATH%\distmat_%LABEL%

@set DIM=1
@set TAUMAX=0

@set NBEGIN=0
@set NEND=0
@set NORM=0
@set SKIPTIME=0

@set DST=%EXP_PATH%\ph_%DATA%_%LABEL%

%BIN%DiffusionRunner_D64.exe --norm %NORM% --nthreads 16 -i %INPUT% -d %DIM% -t %TAUMAX% -o %DST% -b %NBEGIN% -e %NEND% -s %SKIPTIME%
