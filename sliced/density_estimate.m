function im = density_estimate(X,N,sig)
%
% Usage : im = density_estimate(X,N,sig)
%
% Use : Estimate 2D density from point cloud X ([2xn] matrix) :
% 1) X is normalized in [0,1]^2 and then quantified on the [NxN] pixel grid
% 2) im is a N-by-N pixel picture in [0,1] s.t. 
%       im = G_sig * mu_x 
%   where : '*' is the convolution operator
%           'G_sig' is a gaussian kernel with parameter 'sig' and 
%           'mu_x' is the sum of dirac mass prescribed by 'X'
%
% Example 1 : figure, imagesc(density_estimate(rand(2,100)))
% Example 2 : figure, imagesc(density_estimate(randn(2,1e5),2^8,10)), colormap gray, axis equal, axis off
%
% [X,Y] = meshgrid(linspace(0,1,100)); Z = [X(:)';Y(:)']; figure, imagesc(density_estimate(Z,100,10))
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
if nargin<=2
    % define intra distance parameter for density estimation
    if 0 % compute statistics on inner distance
        D = repmat(sum(Z.^2,1),[n 1]); D = D + D' - 2*Z'*Z; % inner-distance
        D = sqrt(D) + diag(inf(n,1)); % valeur 'inf' sur la diagonale
        mini = min(D(:)); % plus petite distance entre deux points
        q10  = quantile(D(D<inf),10); % dist median
        sig = q10*N; % ecart-type à fixer selon la densité de points
    else % just a guess
        sig = 0.02 * N;
    end
end

%% step 1 : normalization in [0,1]^2
Z = ( X - repmat(min(X,[],2), [1 n]) ) ./ repmat( max(X,[],2) - min(X,[],2), [1 n]);

if 0 % add a small border for convolution (5%)
    p = 0.05;
    Z = Z*(1-2*p) + p;
end

% figure, plot(Z(1,:), Z(2,:),'x')

%% step 2 : build kernel

gaussian = @(x,mu,sig) exp(-(x-mu).^2/2/sig^2);
x = 3*round(sig)+1; x = (-x:x);
g = gaussian(x,0,sig);
G = g'*g; G = G/sum(G(:)); % figure, imagesc(G)


%% step 3 : convolution

im = zeros(N); % size(im)
Zi = round(Z*N); % figure, plot(Zi(1,:), Zi(2,:),'x')
di = Z*N - Zi;
Zi(Zi<1) = 1; Zi(Zi>N) = N; 
% figure, plot(Zi(1,:), Zi(2,:),'x')

% interpolation
if 0 % NN interp
    idx = sub2ind([N N], Zi(2,:), Zi(1,:));
    % im(idx) = 1; UNCORRECT if multiple indentical entries
    % use this code instead :
    for i=1:n
        im(idx(i)) = im(idx(i)) + 1;
    end
else % interpolation bi-linéaire
    idxy = sub2ind([N N], Zi(2,:),          Zi(1,:));
    idxY = sub2ind([N N], Zi(2,:),          min(Zi(1,:)+1,N) );
    idXy = sub2ind([N N], min(Zi(2,:)+1,N), Zi(1,:));
    idXY = sub2ind([N N], min(Zi(2,:)+1,N), min(Zi(1,:)+1,N));
    
    for i=1:n
        im(idxy(i)) = im(idxy(i)) +    di(1,i)  .*     di(2,i);   %    xy
        im(idxY(i)) = im(idxY(i)) +    di(1,i)  .*  (1-di(2,i));  %    x(1-y)
        im(idXy(i)) = im(idXy(i)) + (1-di(1,i)) .*     di(2,i);   % (1-x)y
        im(idXY(i)) = im(idXY(i)) + (1-di(1,i)) .*  (1-di(2,i));  % (1-x)(1-y)
    end
end
% figure, imagesc(im)

% convolution
im = conv2(im,G,'same');
% normalization
im = im ./ max(im(:));

return

% TEST
% im0 = 1-im;
% figure, imagesc(im0), colormap gray, axis image, axis off,

