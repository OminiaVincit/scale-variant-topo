@set DATA=BZR

@set BIN=ph-compute\ScaleVariantTopo\Release\
@set EXPNAME=exp_20190516
@set EXPPATH=F:\Research\ScaleVariant\%EXPNAME%\%DATA%

@set KERF=kernel
@rem @set KERF=riemann

@set METHOD=0

@set THRES=0.0
@set INFVAL=0.0
@set TMAX=0.0

@set SCALE=-1
@set T1=0.0
@set T2=1.0

@set PHNAME=ph_20180628_norm_0

@set DIM=1

@set LEFT=%EXPPATH%\barlist_%PHNAME%_d_%DIM%.txt
@set OUT=%EXPPATH%\%KERF%\%PHNAME%_d_%DIM%_method_%METHOD%_T1_%T1%_T2_%T2%_tmax_%TMAX%_thres_%THRES%_inf_%INFVAL%.txt

%BIN%DiagramDistanceRunner_D64.exe --max_tau %TMAX% --dim %DIM% --T1 %T1% --T2 %T2%  -l %LEFT% -o %OUT% --thres %THRES% --method %METHOD% --infval %INFVAL%

@set DIM=0

@set LEFT=%EXPPATH%\barlist_%PHNAME%_d_%DIM%.txt
@set OUT=%EXPPATH%\%KERF%\%PHNAME%_d_%DIM%_method_%METHOD%_T1_%T1%_T2_%T2%_tmax_%TMAX%_thres_%THRES%_inf_%INFVAL%.txt

%BIN%DiagramDistanceRunner_D64.exe --max_tau %TMAX% --dim %DIM% --T1 %T1% --T2 %T2%  -l %LEFT% -o %OUT% --thres %THRES% --method %METHOD% --infval %INFVAL%
