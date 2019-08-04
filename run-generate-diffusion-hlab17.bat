@set MATPATH=F:\tran\Research\ScaleVariant\exp_20190516\REDDIT-BINARY\mat
@set OUTPATH=F:\tran\Research\ScaleVariant\exp_20190516\REDDIT-BINARY\distmat

python generate_diffusion_mat.py --matpath %MATPATH% --outpath %OUTPATH% --normflag 0 --avgflag 0 --taubg 1.0 --taued 50.1 --interval 1.0 --idxbg 1 --idxed 1 --nproc 20