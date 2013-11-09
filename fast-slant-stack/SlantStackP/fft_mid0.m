%  fftmid0 -- 1d fft with argument [-pi,pi] (instead of [0,2pi])
%  Usage:
%    Y = fft_mid0(X)
%  Inputs:
%    X	Array(n) 
%  Outputs:
%    Y   Array(n)
%  Description:
%   Performs 1d fft with grid (-n/2):(n/2-1) on both time
%   and frequency side. 
%     y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
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