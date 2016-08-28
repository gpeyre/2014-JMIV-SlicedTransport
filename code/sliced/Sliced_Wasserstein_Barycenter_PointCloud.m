function [X,E_Y,G] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w,options)

% Sliced Wasserstein Barycenter - multi-dimensional barycenter using sliced wasserstein energy descent
%
%   [X,E_Y,G] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w,options);
%
% OUTPUTs
%
%   X is the Sliced Wasserstein barycenter of Y, minimizing the following energy:
%                   E_Y (X) = Sum_{k} w(k) * SW2(X, Y{k})
%   where SW2(X,Y) = min_{s} 1/2 * Sum_{o,i} { < X[i]-Y[s(i)] , o >^2 } is
%   the quadratic Sliced Wasserstein Distance
%
%   'E_Y' is the energy E_Y(X{k})/N at iterations k = 1 .. niter
%
%   'G' is the average grad norm at each iteration : <||X{k+1} - X{k}||_2>
%   
% INPUTs
%
%  X0 : initialization of barycenter X (should be one of Y distributions)
%
%  Y : list of N-point-wise distributions in R^D (cell array of D-by-N matrices)
%
%  options:
%   options.niter is the number of iterations (100 by default)
%   options.ndir is the number of directions (ndir=dimension by default)
%   options.step is the step descend (step=1 by default).
%   options.eps sets the threshold on normalized displacement/descent amount G, in order to stop prematurely the algorithm
%
%   options.ortho_base=1 uses Gram-Schmidt algorithm for the first 'd' vectors of the base (=dimension) 
%   options.random_base=1 uses random base for each step (mandatory for stochastic gradient descent)
%           else if options.base = [d x ndir] matrix, it is used as fixed Bases
%
%   options.method is the descent strategy : 
%     - 'grad' : gradient descent (classical averaging batch descent which can be combined with Hessian 'normalization')
%     - 'stochastic' : stochastic gradient descent on a subset of directions
%
%   options.hessian=1 uses hessian matrix to scale gradient step
%
%   options.display=1 displays intermediate states (2D/3D only)
%   options.message=1 allows printing
%   options.pause : = 1 adds a pause, = t adds extra time of t seconds at each iteration, 
%
%
% EXAMPLE
%
%     n = 1e3;
%     x = linspace(0,1,n)*2*pi;
%     Y = { [cos(x); sin(x); ]; [cos(x)+4; sin(2*x); ]; [cos(2*x)+2; sin(3*x)+2; ]; };
%     nb = figure; plot(Y{1}(1,:), Y{1}(2,:), 'r'), hold on
%     plot(Y{2}(1,:), Y{2}(2,:), 'g'), plot(Y{3}(1,:), Y{3}(2,:), 'b')
%
%     X0 = Y{1};
%     w = [1 1 1];
%     d = 1e2;
%     t = (1:d)/d*pi;
%     options.ndir = d;
%     options.base = [cos(t); sin(t)];
%     options.method = 'grad';
%     options.ortho_base = 1;
%     options.step = 1;
%     options.hessian = 1;
%     options.niter = 1e2;
%     options.display = 1;
%     tic, [X,E_Y,G] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w,options); toc 
%     figure(nb), plot(X(1,:), X(2,:), '.k')
%
%
% Copyright, J. Rabin & G. Peyré , 2010.
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
nsubdir = getoptions(options, 'nsubdir', ndir); % number of directions
dt = getoptions(options, 'step', 1); % step size
flag_ortho_base = getoptions(options, 'ortho_base', 1); %'d' among 'ndir' directions are orthogonal
Base = getoptions(options, 'base', 0);
flag_method = getoptions(options, 'method', 'grad'); % 'grad' or 'stochastic'
flag_hessian = getoptions(options, 'hessian', 1);

flag_message = getoptions(options, 'message', 1);
flag_display = getoptions(options, 'display', 0);
flag_pause = getoptions(options, 'pause', 0);
epsi = getoptions(options, 'eps', 1e-10);

Vect = @(x) x(:); % simply perform x=X(:); very useful for inline functions


%% Initialization

if dt>1
  disp('/!\ Warning ! iteration step is bigger than one. Continue ?')
  pause
end

X = X0; % initial position
G = []; % descent amount (normalized displacement between two consecutive iterations)
E_Y = []; % Energy at each iteration
flag_end=0; flag_conv=0;

% directions setting (todo : generate minimum correlated directions with blue noise)
if min(size(Base) == [d ndir]) % use the input base
    D = reshape(Base,[d 1 ndir]); clear Base;
else % random basis
    D = randn(d, 1, ndir); % isotropic distribution
end
D = D ./ repmat( sqrt(sum(D.^2,1)), [d 1] ); % normalization to get unit vector

if flag_ortho_base  % (Gram-Schmidt)
    % idem but slower
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

if flag_hessian && strcmpi(flag_method,'grad') % precompute inverse matrix for block descend
    DD = reshape(D(:,1,:),[d ndir]); DD = pinv(DD*DD'); % 'pinv' au lieu de 'inv' permet de traiter les cas ndir<d
end

  

if flag_display % && d<=3
    %random numbers for figures
    fig_num1 = floor(1e3*abs(randn)+1);
    
    color_point = hsv(K);
    figure(fig_num1), hold on,
    for it=1:K
       plot(Y{it}(1,:), Y{it}(2,:), 'o', 'color', color_point(it,:)) 
    end
    
    % bounding box for display
    bb = axis; 
end

%% Precomputation
if strcmpi(flag_method,'grad')
    % sort input points according to selected orientations
    v2{K}=0;
    for kk=1:K % for each point-cloud
        v2{kk} = sum( D.*repmat(Y{kk}, [1 1 ndir]) ); [v2{kk},I2] = sort(v2{kk},2);
    end
    clear I2
end
dx = zeros(d,n,ndir);

%% Algorithm
for it=1:niter 
    if it==niter, flag_end=1; end
    if flag_message, progressbar(it,niter); end
    
     
      
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if strcmpi(flag_method,'grad') % deplacement MOYEN selon le gradient
      
        X_prev = X;

        %%%%% projection + sort (optimal assignment according each direction)
        v1 = sum( D.*repmat(X, [1 1 ndir]) );  [v1,I1] = sort(v1,2);
        % we should use previous sorting 'I1' here

        delta=0; E_Y(end+1)=0;
        
        for kk=1:K % for each point-cloud
            
            if 1 % sortings of input point-clouds are precomputed
                V2 = v2{kk};
            else
                V2 = sum( D.*repmat(Y{kk}, [1 1 ndir]) ); [V2,I2] = sort(V2,2);
            end
            
            for k=1:ndir %pour chaque direction, calcul de l'assignement
                dx(:,I1(1,:,k),k) = repmat(V2(:,:,k)-v1(:,:,k), [d 1]) .* D(:,:,k);        
            end
            delta = delta + w(kk)*sum(dx,3);

            E_Y(end) = E_Y(end) + w(kk) * sum((V2(:)-v1(:)).^2);

        end
        
        E_Y(end) = E_Y(end)/2/ndir;

        if ~flag_hessian,
            dx = delta/ndir; % une normalisation par le nombre de direction est nécessaire
        else
            dx = DD*delta; %avec inverse de la hessienne, pas besoin de diviser par 'ndir'
        end

        X = X + dx*dt;

        G(end+1) = norm(X(:)-X_prev(:),2)/sqrt(n); %deplacement moyen entre deux itérations
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elseif strcmpi(flag_method,'stochastic') % stochastic displacement
       
        X_prev = X;
        
        % random orientation selection
        D = randn(d, 1, ndir); % isotropic distribution
        D = D ./ repmat( sqrt(sum(D.^2,1)), [d 1] ); % normalization to get unit vector

        if flag_ortho_base  % (Gram-Schmidt)
            for k=2:min(d,ndir) %only 'd' orthonormal vectors maximum
                for j=1:k-1
                    D(:,1,k) = D(:,1,k) - D(:,1,j) * (D(:,1,k)' * D(:,1,j));
                end
                D(:,1,k) = D(:,1,k)/sqrt(sum(D(:,1,k).^2,1));
            end
        end
        D = repmat(D,[1 n 1]);
        

        %%%%% projection + sort (optimal assignment according each direction)
        v1 = sum( D.*repmat(X, [1 1 ndir]) );  [v1,I1] = sort(v1,2);
        
        delta=0; E_Y(end+1)=0;
        
        for kk=1:K %pour chaque distribution
            
            % sort values of Y's
            V2 = sum( D.*repmat(Y{kk}, [1 1 ndir]) ); [V2,I2] = sort(V2,2);
            
            % dx = zeros(d,n,ndir); % useless
            for k=1:ndir %pour chaque direction, calcul de l'assignement
                dx(:,I1(1,:,k),k) = repmat(V2(:,:,k)-v1(:,:,k), [d 1]) .* D(:,:,k);        
            end
            delta = delta + w(kk)*sum(dx,3);

            E_Y(end) = E_Y(end) + w(kk)*norm( V2(:) - v1(:),2).^2;

        end

        if ~flag_hessian,
            dx = delta/ndir; % une normalisation par le nombre de direction est nécessaire
        else % using hessian normalization
            DD = reshape(D(:,1,:),[d ndir]); DD = pinv(DD*DD'); % /!\ [dxd] inverse matrix computed at each step !
            dx = DD*delta;
        end

        X = X + dx*dt;

        G(end+1) = norm(X(:)-X_prev(:),2)/sqrt(n); %deplacement moyen entre deux itérations
    
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % Display energy decrease and point's moves
    if flag_display && (mod(it*10,niter)==0 || it==1)
        figure(fig_num1); clf;
        
        subplot(1,2,1), hold on;
        for ii=1:K
           plot(Y{ii}(1,:), Y{ii}(2,:), 'o', 'color', color_point(ii,:)) 
        end
           plot(X(1,:), X(2,:), '.k') 
        title(['Iteration ' num2str(it) '/' num2str(niter)]);
        hold off; %axis equal
        axis(bb);
        %drawnow;
        
        
        %error plot display
        figure(fig_num1); subplot(1,2,2)
        semilogy(E_Y/E_Y(1),'b')
        hold on;
        semilogy(G/G(1),'r')
        %legend('Energy','Norm') = ralentissement !
        if epsi>0, axis([0 niter epsi/G(1) 1]); 
        else axis([0 niter eps/G(1) 1]); end %min(G(:))
        if epsi>0, line([0 1000],[epsi epsi]/G(1),'color','m'); end
        drawnow;
    end
    
    if flag_pause==1, pause, elseif flag_pause, pause(flag_pause), end
    
    if G(it)<=epsi
        if flag_message, 
            disp(' '), disp(['= convergence reached (' num2str(it) ' iterations)']); 
            disp(['= Barycenter Energy : ' num2str(E_Y(end)) ]);
        end
        
        flag_end=1; flag_conv=1;
    end
    
    if flag_end, break, end
end

if flag_display 
    figure(fig_num1); subplot(1,2,2)
    legend('SW2 energy','gradient norm',4)
end

if ~flag_conv
    if flag_message, 
      disp(' '),
      disp(['> convergence not reached after ' num2str(it) ' iterations']);
      disp(['> Last move (average grad norm : <||X{k+1} - X{k}||_2>) : ' num2str(G(end)) ]);
      disp(['> Barycenter Energy : ' num2str(E_Y(end)) ]);
      disp(' '),
    end
end    



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

