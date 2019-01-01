@set BIN=ph-compute\ScaleVariantTopo\Release\

@set TIMEHOLE=0.0
@set TIMETAU=1.0

@set THRES=0.0
@set METHOD=0
@set TAUMAX=0

@set LEFT=demo\GN_barcodes_dim_%DIM%.txt
@set OUT=demo\kernel\GN_dim_%DIM%_kernel_method_%METHOD%_timehole_%TIMEHOLE%_timetau_%TIMETAU%_thres_%THRES%_taumax_%TAUMAX%.txt

%BIN%DiagramDistanceRunner_D64.exe --dim %DIM% --timehole %TIMEHOLE% --timetau %TIMETAU%  -l %LEFT% -o %OUT% --thres %THRES% --method %METHOD% --taumax %TAUMAX%

