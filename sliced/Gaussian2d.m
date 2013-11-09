function mat = Gaussian2d(...
    gsize, center, sigmax, sigmay, theta, offset, factor)
% GAUSSIAN2d - Produces a 2-dimensional gaussian matrix
%
%   Input
%   gsize   :  Size of the output 'gauss', should be a 1x2 vector
%   center  :  The center position of the gaussian [row col]
%   sigmax  :  Std. dev. in the X direction
%   sigmay  :  Std. dev. in the Y direction
%   theta   :  Rotation in degrees
%   offset  :  Minimum value in output
%   factor  :  Related to maximum value of output, should be 
%                    different from zero  
%
%   Output
%   mat     : probability vector for each angle bin 
%
%   Author : Nikolas Markou 
%
%   Exemple : n = 100; G = Gaussian2d([n n], [n n]/2, 20, 10, 45, 0, 1); figure, imagesc(G)
%-------------------------------------------------------------------------
%%  Validate input arguments
    if ~exist('gsize','var')
        error('Error : gsize argument not set');
    end
    
    if ~exist('center','var')
        center = gsize / 2;
    end
    
    if ~exist('sigmax','var')
        sigmax = 1;
    end
    
    if ~exist('sigmay','var')
        sigmay = 1;
    end
    
    if ~exist('theta','var')
        theta = 0;
    end
     
    if ~exist('offset','var')
        offset = 0;
    end
    
    if ~exist('factor','var')
        factor = 1;
    end
    theta  = (theta/180)*pi;  
%-------------------------------------------------------------------------
%%  Compute
    [X,Y] = meshgrid(1:gsize(1),1:gsize(2));
    Xc  = X - center(1);
    Yc  = Y - center(2);
    xm  = (Xc)*cos(theta) - (Yc)*sin(theta);
    ym  = (Xc)*sin(theta) + (Yc)*cos(theta);
    u   = (xm / sigmax).^2 + (ym / sigmay).^2;
    mat = offset + factor * exp(-u/2);
end
