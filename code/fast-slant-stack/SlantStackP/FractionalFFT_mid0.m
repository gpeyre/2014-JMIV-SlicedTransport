%  FractionalFFT_mid0: Fractional Fourier Transform (midpoint at 0)
%   Usage:
%     xhat = FractionalFFT_mid0(x,alpha)
%   Inputs:
%      x       array(n), n dyadic
%      alpha   scalar -- fractional exponent   
%   Outputs:
%      xhat    array(n), n dyadic
%   Description:
%      Performs fft of series x with an nonstandard "alpha" factor
%   in exponent, so matrix element W_{j,k} = e^(-alpha*i*2*pi*k*t/n).
%       alpha=1   corresponds to ordinary fft.
%       alpha=-1  corresponds to (unnormalized) ordinary ifft.
%   Transforms with +/-alpha are adjoints of each other
%   Both the inputs and outputs are centered so that
%   k and t run through -n/2 <= k,t < n/2 
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