%  Pseudopolar_double: Pseudopolar FFT
%   Usage:
%      Y=Pseudopolar_double(X)
%   Inputs:
%     X      n*n matrix (x,y)
%   Outputs:
%     Y      2n*2n matrix (theta,r)
%   Description:
%     Performs pseudopolar
%   FFT on n*n matrix by using frft.
% 
%  Performs pseudopolar FFT on square n*n matrix by using frft.
%  The resulting 2n*2*n matrix is divided into quadrants
%   -----
%   -----           -------
%  | 1 /||        | 1 | 3 |  Y(1,:) is  diagonal line 
%  |4 /2||  - Y= |___|___|  Y(:,1) is _  angles
%  | /  ||  -/    | 2 | 4 |             |
%  |/ 3 ||        |   |   |
%   -----           -------
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