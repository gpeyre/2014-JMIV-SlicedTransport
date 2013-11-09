%  UFWT_CDJV -- Forward Wavelet Transform (boundary-corrected, NOT preconditioned)
%   Usage
%     wc = UFWT_CDJV(x,L,N)
%   Inputs
%     y    1-d signal, length(x) = 2^J
%     L    Level of V_0;  L << J
%     N    Degree of Daubechies Filters
% 
%   Description
%     CDJV have developed an algorithm for wavelets on the interval which
%     preserves the orthogonality, vanishing moments, smoothness, and compact
%     support of Daubechies wavelets on the line.
% 
%     The algorithm for wavelets on the interval of CDJV involves four objects
%     not present in the usual periodized algorithm: right edge filters, left
%     edge filters, and pre- and post- conditioning operators.
% 
%     These objects are supplied by appropriate requests to MakeCDJVFilter.
% 
% 	 This variant does not apply the preconditioning operator.
% 
%     To reconstruct use CDJV_IWT.
% 
%   See Also
%     IWT_CDJV, FWT_PO, IWT_PO, MakeCDJVFilter
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