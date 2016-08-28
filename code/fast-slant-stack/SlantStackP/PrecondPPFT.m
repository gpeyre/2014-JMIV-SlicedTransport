%  PrecondPPFT: Preconditioner for Pseudopolar FFT
%   Usage:
%    PrPPFT = PrecondPPFT(PPFT,Pre,Pow)
%   Inputs:
%    PPFT		array(m,m) PseudopolarFFT
%    Pre			2/1 Preconditioner type
%    Pow         exponent of preconditioner e.g. 1,1/2,-1/2,-1
%   Outputs:
%    PrPPFT      array(m,m) Preconditioned PseudopolarFFT
%   Description:
%    The Pseudopolar FFT array is multiplied by a function of
%     pseudo-radius which is a chosen power of pseudo radius.
%   See Also
%    Pseudopolar FFT, Adj_Pseudopolar FFT
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