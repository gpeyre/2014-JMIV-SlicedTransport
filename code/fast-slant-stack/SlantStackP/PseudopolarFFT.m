%  PseudopolarFFT: Pseudopolar FFT
%   Usage:
%      Y = PseudopolarFFT(X,Pre)
%   Inputs:
%     X      n*n matrix (x,y)
%     Pre	Preconditioning flag (2/1/0) default=0
%   Outputs:
%     Y      2n*2n matrix (theta,r)
%   Description:
%     Performs pseudopolar FFT on n*n matrix.
% 
%  Performs pseudopolar FFT on square n*n matrix..
%  The resulting 2n*2*n matrix has two `panels'
%  	An n * 2n panel for `mostly horizontal lines'...
%    An n * 2n panel for `mostly vertical lines'...
%  We have the following association:
%  Row			Slopes of Level Sets of sinusoids
%  Y(1,:)         mostly horiz, slope +1  
%  Y(n/2+1,:)     horiz 
%  Y(n,:)         mostly horiz, slope -(n-1)/n
%  Y(n+1,:)       mostly vert, slope -1
%  Y(3n/2+1,:)    vert
%  Y(2n,:)        mostly vert, slope (n-1)/n
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