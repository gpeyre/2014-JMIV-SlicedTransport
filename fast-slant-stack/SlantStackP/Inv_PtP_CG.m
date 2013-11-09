%  Inv_PtP_CG: Generalized inverse of Pseudopolar FFT Frame
%   Usage:
%      Y = Inv_PtP_CG(X,Age,Pre,MaxIts,ErrTol)
%   Inputs:
%     X      n*n matrix (theta,r)
%     Age	Flag (1/0) for use of old code default=0
%     Pre	Preconditioning flag (2/1/0) default=1
%     MaxIts max number of iterations, default 5
%     ErrTol error tolerance, default 1.e-9
%   Outputs:
%     Y      n*n matrix (x,y)
%   Description:
%     Performs approximate inverse of P'P using conjugate
%     gradients. Usually requires 3 or 4 steps only.
%   See Also
% 	PtP, PseudopolarFFT, Adj_PseudopolarFFT
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