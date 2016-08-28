%  UIWT_CDJV -- Inverse Wavelet Transform  (boundary corrected)
%   Usage
%     x = UIWT_CDJV(wc,L,N)
%   Inputs
%     wc   1-d wavelet transform
%     L    Level of V_0;  L << J
%     N    Degree of Daubechies Filters
%   Outputs
%     x    1-d signal: length(y) = 2^J
% 
%   See Also
%     FWT_CDJV, MakeCDJVFilter
% 
%   References
%    This is an implementation of the Cohen-Daubechies-Jawerth-Vial Algorithm
%    for orthonormal wavelet bases of compact support, with boundary corrected
%    wavelets at 0 and 1.
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