%% script to compare the sliced wasserstein formulation for point clouds with exact wasserstein projection
% julien.rabin@unicaen.fr - 2013

clear all, 
%close all,
%clc

% parameters
d = 2; % number of directions
n = 1e2;

Nmax_comp = 2e3; % maximum number of points for comparison with Mosek
Nmax_plot = 5e3; % maximum number of points for color point display


%% 0) Generate input point-clouds

x = linspace(0,1,n)*2*pi;
Y = [];
if 1 % generate curves

    % Y{end+1} = [cos(4*x); sin(4*x) ];  % circle
    % Y{end+1} = [cos(x)+.1*cos(15*x); sin(x)+.1*sin(15*x) ];  % flower
    % Y{end+1} = [cos(x)+.2*cos(15*x); sin(x)+.2*sin(15*x) ];  % flower
    % Y{end+1} = [cos(x); .5*sin(x); ]; % ellipse
    % Y{end+1} = [.5*cos(x)+2; sin(x); ]; % ellipse
    % Y{end+1} = [ 1.2*sin(2*x); .5*cos(x); ]; % 8
    % Y{end+1} = [.5*cos(x); 1.2*sin(2*x); ]; % infinity
    % Y{end+1} = [cos(2*x).^3+4; sin(3*x).^3; ]; % arrow
     Y{end+1} = [cos(x).^3; sin(x).^3 ];  % astroid
     t = linspace(-1,1,n/2); y = t.^2; Y{end+1} = - [ (y+1)/2+0, 1-sqrt(1-y)+0; t,t ]; % Crescent

    % symmetry
    % Y{end+1} = - Y{end};
     
elseif 1

    N = n*3*4*5; % required to get polygone figures
    X = {};

    a = 108/180*pi; % inner pentagone angle
    b = pi-a; cb = cos(b); sb = sin(b);
    c = b - (pi/2 - b); cc = cos(c); sc = sin(c);
    L = 5;

    % circle
    t = (1:N)/N*2*pi;
    X{1} = [cos(t); sin(t)]; % in [-1,1]^2
    
    % square
    n = N/4;
    t = (0:n-1)/n;
    z = zeros(1,n);
    o = z+1;
    X{2} = [ [t;z] , [o;t] , [1-t;o] , [z;1-t] ]*2-1; % in [-1,1]^2
    
    % equilateral triangle
    n = N/3;
    t = (0:n-1)/n;
    z = zeros(1,n);

    X{3} = [ [t;z] , [(2-t)/2 ; t*sqrt(3)/2] , [(1-t)/2 ; (1-t)*sqrt(3)/2] ];
    X{3} = X{3}*2+repmat([-1;-sqrt(3)/2],[1,N]); % in [-1,1]^2
    
    % pentagone
    n = N/5;
    t = (0:n-1)/n;
    z = zeros(1,n);

    X{4} = [ [(t-1)*cb;(1-t)*sb] , [t;z] , [1+t*cb ; t*sb] , [1+cb - t*sc ; sb + cc*t] , [1+cb-sc - t*sc ; sb + cc *(1-t)] ];
    X{4} = X{4} + repmat([-1/2;-sb],[1,N]); % in [-1,1]^2
    
    % star
    o = z+1;
    X{5} = [ [t*(1+cb);t*sb] , [1+cb-t*(1+2*cb);o*sb] , [t*(1+cb)-cb;(1-t)*sb] , [1 + t*(cb-sc) ; (sb + cc)*t] , [(1+cb-sc) * (1 - t); (sb + cc) *(1-t)] ];
    X{5} = X{5} + repmat([-1/2;-sb],[1,N]); % in [-1,1]^2
    %X{5}(1,:) = X{5}(1,:) + (1+cb-sc);
    %X{5}(2,:) = X{5}(2,:) + (sb+cc);

    % select two figures
    Y = {X{4}; X{5}};
    n = N; % Number of points
    
else % mixture of two gaussians
    n1 = round(n/3);
    n2 = n-n1;

    theta = [1/3,0; 2/3,1/4; 1/2,1]*pi % (0:K-1)/K*pi + pi/3; % angle
    sx = [1,1; 1,.5; ]/4; % first  standard deviation
    sy = [1,1; .5,1; ]/4; % second standard deviation
    mux = [-2,0; 1,-1];
    muy = [-2,0; 0,0];
    
    for j = 1:2
        x1 = randn(1,n1) * sx(j,1); x2 = randn(1,n2) * sx(j,2);
        y1 = randn(1,n1) * sy(j,1); y2 = randn(1,n2) * sy(j,2);
        t1 = theta(j,1); t2 = theta(j,2);
        c1 = cos(t1); s1 = sin(t1); c2 = cos(t2); s2 = sin(t2);
        Y{j} = [ c1*x1-s1*y1 + mux(j,1), c2*x2-s2*y2 + mux(j,2); ...
                 s1*x1+c1*y1 + muy(j,1), s2*x2+c2*y2 + muy(j,2); ...
               ];
    end
end

% translation ou scaling
if 1
    Y{2}(1,:) = Y{2}(1,:) + 4;
    Y{2}(2,:) = Y{2}(2,:) + 0;
else
    Y{1} = Y{1}*2;
    Y{2} = Y{2}/2; 
end

if 0 % TEST use random permutation
    for i=1:numel(Y)
        Y{i} = Y{i}(:,randperm(n));
    end
end

nb = figure; 
if n > Nmax_plot
    plot(Y{1}(1,:), Y{1}(2,:), 'or'), hold on
    plot(Y{2}(1,:), Y{2}(2,:), 'og')
else
    if (n <= Nmax_comp)
        subplot(1,2,1)
    end
    colour = jet(n);
    plot(Y{2}(1,:), Y{2}(2,:), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize',10),
    hold on
    D = Y{1} - repmat(mean(Y{1},2), [1 n]);
    [D,i] = sort(sum(D.^2));
    Y{1} = Y{1}(:,i);
    scatter(Y{1}(1,:), Y{1}(2,:),30,colour,'fill')
end
axis equal  
  
%% 1) Compute SW2 projection from 2 point-clouds
  
% parameters setting
X0 = Y{1};
w = [0 1]; w = w/sum(w); % weights for projection
t = (1:d)/d*pi; % uniform sampling or orientations
options.ndir = d;
options.base = [cos(t); sin(t)]; % directions set
options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'stochastic'; % 'grad' or 'stochastic'
options.ortho_base = 1;
options.niter = 1e3;
options.epsi = 1e-3;
options.display = 1;

tic, [SW2_proj,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w,options); toc

figure(nb),
if n > 1e3
    plot(SW2_proj(1,:), SW2_proj(2,:), '.b')
else
    scatter(SW2_proj(1,:), SW2_proj(2,:),30,colour,'fill')
    title(['SW2 projection with ' num2str(options.ndir) ' directions'])
    legend('Y1', 'Y2', 'SW2_proj')
end
drawnow

W2_approx = sum((SW2_proj(:)-X0(:)).^2) % display approx L2-OT cost
      
%% 2) comparison with exact OT (bipartite graph matching with Hungarian Algorithm or Mosek library)

if (n>Nmax_comp)
    disp('STOP : computation may be very long with such number of points !')
    
    return
else
    % build cost matrix : Cij = || Y1(i) - Y2(j) ||^2
    Y12 = repmat(sum(Y{1}.^2,1) , [n 1]); % L2^2 norm of Y1
    Y22 = repmat(sum(Y{2}.^2,1)', [1 n]); 
    costMat = Y12 + Y22 - 2*Y{2}'*Y{1}; % L2^2 distance : figure, imagesc(costMat)
    
    if (n<=2e2) % Hungarian Algorithm
        disp('Hugarian algorithm!')
        tic, 
        [perm,W2] = munkres(costMat');
        toc

    elseif (n<=Nmax_comp) % Mosek
        disp('warning : computation may be very too long with Hugarian algorithm!')
        disp('we use Mosek library instead (faster)')

        % to get a free trial licence, go to : http://license.mosek.com/license2/academic/
        
        [W2,F,perm] = OT_mosek(costMat');
    end
end

W2 % display exact L2-OT cost
W2_proj = Y{2}(:,perm);

figure(nb),
if n > Nmax_plot
    plot(W2_proj(1,:), W2_proj(2,:), 'om')
    title('Wasserstein projection and its sliced approximation')
    legend('Y1', 'Y2', ['SW2_proj, ' num2str(options.ndir) ' dir'], 'W2_proj')
else
    subplot(1,2,2)
    plot(Y{2}(1,:), Y{2}(2,:), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize',10),
    hold on
    scatter(Y{1}(1,:), Y{1}(2,:),30,colour,'fill')
    
    scatter(W2_proj(1,:), W2_proj(2,:),30,colour,'fill')
    title('W2 projection')
    legend('Y1', 'Y2', 'W2_proj')
    
    axis equal  
end
drawnow


%% 3) Compare point-clouds barycenters

rho = .5; % in [0,1]

 W2_bary = (1-rho)*X0 + rho* W2_proj;
SW2_bary = (1-rho)*X0 + rho*SW2_proj;

% compute the weighted barycenter by interpolating the matching flow
figure, 
plot(Y{1}(1,:), Y{1}(2,:), '.b'), hold on
plot(Y{2}(1,:), Y{2}(2,:), '.k')
plot(SW2_bary(1,:),SW2_bary(2,:), 'xg')
plot( W2_bary(1,:), W2_bary(2,:), 'om')
title('Wasserstein interpolation')
legend('Y1', 'Y2', ['SW2 interp with' num2str(options.ndir) ' dir'], 'W2 interp')

%% 4) density estimation and comparison

% convolution with gaussian kernel
N = 2^8;
d_SW2_bary = 1 - density_estimate(SW2_bary,N,N/100*4);
d_W2_bary  = 1 - density_estimate( W2_bary,N,N/100*4);

figure, 
subplot(1,2,1), imagesc(d_SW2_bary), title('SW2 interpolation')
colormap gray, axis image, axis off,
subplot(1,2,2), imagesc(d_W2_bary ), title(' W2 interpolation') 
colormap gray, axis image, axis off,
  
if 0 % display difference
    figure,
    imagesc(d_SW2_bary - d_W2_bary)
end
