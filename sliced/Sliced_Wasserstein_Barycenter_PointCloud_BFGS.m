function [X] = Sliced_Wasserstein_Barycenter_PointCloud_BFGS(X0,Y,w,options)

% Sliced Wasserstein Barycenter - multi-dimensional barycenter using sliced
% wasserstein energy descent with BFGS (quasi-Newton minimization algorithm)
%
%   [X] = Sliced_Wasserstein_Barycenter_PointCloud_BFGS(X0,Y,w,options);
%
% OUTPUTs
%
%   X [dxN] is the Sliced Wasserstein barycenter of Y, minimizing the following energy:
%                   E_Y (X) = Sum_{k} w(k) * SW2(X, Y{k})
%   where SW2(X,Y) = min_{s} 1/2 * Sum_{o,i} { < X[i]-Y[s(i)] , o >^2 } is
%   the quadratic Sliced Wasserstein Distance
%
% INPUTs
%
%  X0 : initialization of barycenter X (should be one of Y distributions)
%
%  Y : list of N-point-wise distributions in R^D (cell array of D-by-N matrices)
%
%  options:
%   options.niter is the maximum number of iterations (100 by default)
%   options.ndir is the number of directions (ndir=dimension by default)
%   options.ortho_base = 1 uses Gram-Schmidt algorithm for the first 'd' vectors of the base (=dimension) 
%   if options.base = [d x ndir] matrix, it is used as fixed Bases
%
%
% EXAMPLE
%
%     n = 1e3;
%     x = linspace(0,1,n)*2*pi;
%     Y = { [cos(x); sin(x); ]; [1/2*cos(x) .* (1-cos(x))+4; 1/2*sin(x) .* (1-cos(x)); ]; [cos(x).^3+2; sin(x).^3+2; ]; };
%     nb = figure; plot(Y{1}(1,:), Y{1}(2,:), 'r'), hold on
%     plot(Y{2}(1,:), Y{2}(2,:), 'g'), plot(Y{3}(1,:), Y{3}(2,:), 'b')
%
%     X0 = Y{1};
%     w = [1 1 1];
%     d = 1e2;
%     t = (1:d)/d*pi;
%     options = [];
%     options.ndir = d;
%     options.base = [cos(t); sin(t)];
%     options.niter = 1e2;
%     tic, [X1] = Sliced_Wasserstein_Barycenter_PointCloud_BFGS(X0,Y,w,options); toc 
%     figure(nb), plot(X1(1,:), X1(2,:), '.m')
%     tic, [X2] = Sliced_Wasserstein_Barycenter_PointCloud     (X0,Y,w,options); toc 
%     figure(nb), plot(X2(1,:), X2(2,:), '.k')
%     legend('Y1', 'Y2', 'Y3', 'SW2 quasi-newton', 'SW2 newton')
%     axis equal, 
%
% Copyright, J. Rabin & G. Peyré , 2013.
% CNRS - Cérémade, Dauphine University, Paris, France
%
% See Also SLICED_WASSERSTEIN_PROJECTION, SLICED_WASSERSTEIN_KMEANS
%

%% extract options

n = size(X0,2);
d = size(X0,1);

K = length(Y);
if ~exist('w','var') || length(w)<K
  w=ones(1,K)/K;
else
  w=w/sum(w);
  w = reshape(w,1,K);
end

options.null = 0;
niter = getoptions(options, 'niter', 100);
ndir = getoptions(options, 'ndir', d); % number of directions
flag_ortho_base = getoptions(options, 'ortho_base', 1); %'d' among 'ndir' directions are orthogonal
Base = getoptions(options, 'base', 0);

Vect = @(x) x(:); % simply perform x=X(:); very useful for inline functions

%% Initialization

% directions setting 
if min(size(Base) == [d ndir]) % use the input base
    D = reshape(Base,[d 1 ndir]); clear Base;
else % random basis
    D = randn(d, 1, ndir); % isotropic distribution
end
D = D ./ repmat( sqrt(sum(D.^2,1)), [d 1] ); % normalization to get unit vector

if flag_ortho_base  % (Gram-Schmidt)
    % idem but longer
    %Q = orth(reshape(D(:,1,1:min(d,ndir)),[d min(d,ndir)]));
    %D(:,1,1:min(d,ndir)) = reshape(Q,[d 1 min(d,ndir)]);

    for k=2:min(d,ndir) %only 'd' orthonormal vectors maximum
        for j=1:k-1
            D(:,1,k) = D(:,1,k) - D(:,1,j) * (D(:,1,k)' * D(:,1,j));
        end
        D(:,1,k) = D(:,1,k)/sqrt(sum(D(:,1,k).^2,1));
    end
end

D = repmat(D,[1 n 1]);

%% Precomputation

% sort input points according to selected orientations
v2{K}=0;
for kk=1:K % for each point-cloud
    v2{kk} = sum( D.*repmat(Y{kk}, [1 1 ndir]) ); [v2{kk},~] = sort(v2{kk},2);
end

%% BFGS

% define gradient function
Grad = @(X) gradientSW2(X,Y,w,D,v2);
% use : [E,g] = Grad(f); where E is [1-by-1] and g is [N-by-1]

opts = [];
opts.bfgs_memory = 100;
opts.niter = options.niter;

% BFGS
[X, ~, ~] = perform_bfgs(Grad, X0(:), opts);

X = reshape(X,[d n]);

function [SW2,gradSW2] = gradientSW2(u,Y,w,D,v2)

[d,n] = size(Y{1});
K = numel(Y);
ndir = size(D,3);

X = reshape(u,[d,n]);

%%%%% projection + sort (optimal assignment according each direction)
v1 = sum( D.*repmat(X, [1 1 ndir]) );  [v1,I1] = sort(v1,2);

gradSW2=0; SW2=0;

dx = zeros(d,n,ndir); 
for kk=1:K % for each point-cloud

    V2 = v2{kk};

    for k=1:ndir %pour chaque direction, calcul de l'assignement
        dx(:,I1(1,:,k),k) = repmat(v1(:,:,k)-V2(:,:,k), [d 1]) .* D(:,:,k);        
    end
    gradSW2 = gradSW2 + w(kk)*sum(dx,3);

    SW2 = SW2 + w(kk) * sum((v1(:)-V2(:)).^2)/2;
end

gradSW2 = gradSW2(:);

return;


function v = getoptions(options, name, v, mendatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end 

