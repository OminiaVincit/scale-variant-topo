@set MATPATH=F:\Research\ScaleVariant\exp_20190516\DD\mat
@set OUTPATH=F:\Research\ScaleVariant\exp_20190516\DD\distmat

python generate_diffusion_mat.py --matpath %MATPATH% --outpath %OUTPATH% --normflag 0 --avgflag 0 --taubg 1.0 --taued 100.0 --interval 2.0 --idxbg 1 --idxed 1