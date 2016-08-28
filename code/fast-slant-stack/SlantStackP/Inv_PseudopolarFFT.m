%  Inv_PseudopolarFFT: Generalized inverse of Pseudopolar FFT
%   Usage:
%      X = Inv_PseudopolarFFT(Y,Age,Pre,DC,K,E)
%   Inputs:
%     Y      2n*2n matrix (theta,r)
%     Age	Flag (1/0) for use of old code default=0
%     Pre	Preconditioning flag (1/0) default=0
%     MaxIts max number of iterations. Default 5.
%     ErrTol error tolerance. Default 1.e-9.
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Performs inverse of pseudopolar
%   FFT on 2n*2n matrix by using frft.
%   This is an approximate inverse.
% 
%  This function finds the inverse polar fft
%  by using a conjugate gradient solver on the associated Gram system.
%  The number of iterations is bounded by MaxIts
%  and the residual error by ErrTol.
% 
%
%
% Part of BeamLab Version:200
% Built:Friday,23-Aug-2002 00:00:00
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail beamlab@stat.stanford.edu
%
%
% Part of BeamLab Version:200
% Built:Saturday,14-Sep-2002 00:00:00
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail beamlab@stat.stanford.edu
%