function [im, bbox] = density_estimate2(X,N,sig,bbox)
%
% Usage : im = density_estimate(X,N,sig,box)
%
% Use : Estimate 2D density from point cloud X ([2xn] matrix) :
% 1) 'X' is normalized in 'bbox' [2x2] and then quantified on the [NxN] pixel grid
%     bbox = [xmin, xmax; ymin, ymax] by default
% 2) 'im' is a N-by-N pixel picture in [0,1] s.t. 
%       im = G_sig * mu_x 
%   where : '*' is the convolution operator
%           'G_sig' is a gaussian kernel with parameter 'sig' and 
%           'mu_x' is the sum of dirac mass prescribed by 'X'
%
% Example 1 : figure, imagesc(density_estimate2(rand(2,100)))
% Example 2 : figure, imagesc(density_estimate2(randn(2,1e5),2^8,10,[-1,1;-1,1]*2)), colormap gray, axis equal, axis off
%
% [X,Y] = meshgrid(linspace(0,1,10)); Z = [X(:)';Y(:)']; figure, imagesc(density_estimate2(Z,100))
%
% julien.rabin@unicaen.fr - 2013

[d,n] = size(X);
if d~=2
    disp('X should be [2xn] matrix')
    return
end

% default parameter
if nargin==1
    N = 2^8; % pixel resolution 
end

if numel(N) == 2
    M = N(2); %: number of lines
    N = N(1); % columns
elseif numel(N)==1
    M = N; % default
else
    error('N should be either a scalar or an array with two values !')
end
    
if nargin<=2
    % define intra distance parameter for density estimation
    if 0 % compute statistics on inner distance
        D = repmat(sum(Z.^2,1),[n 1]); D = D + D' - 2*Z'*Z; % inner-distance
        D = sqrt(D) + diag(inf(n,1)); % valeur 'inf' sur la diagonale
        mini = min(D(:)); % plus petite distance entre deux points
        q10  = quantile(D(D<inf),10); % dist median
        sig = q10*sqrt(N*M); % ecart-type à fixer selon la densité de points
    else % just a guess
        sig = 6 * sqrt(N*M) / n;
    end
end

if nargin<=3
    bbox = [ min(X,[],2), max(X,[],2) ]; % [xmin, xmax; ymin, ymax]
end

%% step 1 : normalization in [0,1]^2
Z = ( X - repmat(bbox(:,1), [1 n]) ) ./ repmat( bbox(:,2) - bbox(:,1) , [1 n]);

% figure, plot(Z(1,:), Z(2,:),'x')

%% step 2 : build normalized kernel

gaussian = @(x,mu,sig) exp(-(x-mu).^2/2/sig^2);
x = 3*round(sig)+1; x = (-x:x);
g = gaussian(x,0,sig);
G = g'*g; G = G/sum(G(:)); % figure, imagesc(G)


%% step 3 : convolution

im = zeros(N,M); % size(im)
Zq = Z.*repmat([M;N]-1,[1 n]) + 1;
Zi = floor( Zq );
di = Zq - Zi; % figure, plot(di(1,:), di(2,:),'x')
Zi(Zi<1) = NaN; Zi(1,Zi(1,:) > M) = NaN;  Zi(2,Zi(2,:) > N) = NaN; 
% figure, plot(Zi(1,:), Zi(2,:),'x')

% interpolation
if 0 % NN interp
    idx = sub2ind([N M], Zi(2,:), Zi(1,:));
    % im(idx) = 1; UNCORRECT if multiple indentical entries
    % use this code instead :
    j = (1:n);
    for i = j(~isnan(idx))
        im(idx(i)) = im(idx(i)) + 1;
    end
else % interpolation bi-linéaire
    idxy = sub2ind([N M], Zi(2,:),          Zi(1,:));
    idxY = sub2ind([N M], Zi(2,:),          min(Zi(1,:)+1,M) );
    idXy = sub2ind([N M], min(Zi(2,:)+1,N), Zi(1,:));
    idXY = sub2ind([N M], min(Zi(2,:)+1,N), min(Zi(1,:)+1,M));
    j = (1:n);
    for i = j(~isnan(idxy))
        im(idxy(i)) = im(idxy(i)) +    di(1,i)  .*     di(2,i);   %    xy
        im(idxY(i)) = im(idxY(i)) +    di(1,i)  .*  (1-di(2,i));  %    x(1-y)
        im(idXy(i)) = im(idXy(i)) + (1-di(1,i)) .*     di(2,i);   % (1-x)y
        im(idXY(i)) = im(idXY(i)) + (1-di(1,i)) .*  (1-di(2,i));  % (1-x)(1-y)
    end
end
% figure, imagesc(im)

% convolution
im = conv2(im,G,'same');
% normalization by the number of points
im = im ./ n;
% im = im ./ max(im(:));

return

