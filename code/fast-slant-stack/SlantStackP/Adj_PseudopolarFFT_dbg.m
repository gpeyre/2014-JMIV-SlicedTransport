%  Adj_PseudopolarFFT: Adjoint of Pseudopolar FFT
%   Usage:
%      X = Adj_PseudopolarFFT(Y)
%   Inputs:
%     Y      2n*2n matrix (theta,r)
%     Pre	Preconditioning flag (1/0) default=0
%     DC     Flag for Special Treatment of DC term (1/0) default=0
%   Outputs:
%     X      n*n matrix (x,y)
%   Description:
%     Performs adjoint of pseudopolar FFT.
%   This is an approximate inverse to PseudopolarFFT.
%   See Also
%     PseudopolarFFT, Inv_PseudopolarFFT
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