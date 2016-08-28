%% script to compare the sliced wasserstein formulation for point clouds with exact wasserstein projection
% julien.rabin@unicaen.fr - 2013

addpath('toolbox/');

% clear all,
%close all,
% clc

RANDOM_INIT = 1
RANDOM_DIR = 0
INCREASE_DIR = 1

% parameters
n = 100;
Nmax_comp = 2e3; % maximum number of points for comparison with Mosek


%% 0) Generate input point-clouds

x = linspace(0,1,n)*2*pi;
Y = [];
if 1 % generate curves

    % Y{end+1} = [cos(x)+.1*cos(15*x); sin(x)+.1*sin(15*x) ];  % flower
     Y{end+1} = [cos(x); .5*sin(x); ]; % ellipse
     Y{end+1} = [.5*cos(x)+2; sin(x); ]; % ellipse
    % Y{end+1} = [ 1.2*sin(2*x); .5*cos(x); ]; % 8
    % Y{end+1} = [.5*cos(x); 1.2*sin(2*x); ]; % infinity
    % Y{end+1} = [cos(2*x).^3+4; sin(3*x).^3; ]; % arrow
    
    
    a = 45/180*pi; ca = cos(a); sa = sin(a);
    Y{end} = [ca,-sa;sa,ca] * Y{end};

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
    a = 45/180*pi; ca = cos(a); sa = sin(a);
    Y = {[ca,-sa;sa,ca] * X{3}; X{4}};
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

% translation
Y{2}(1,:) = Y{2}(1,:) + 6;
clear X
  

if 0 % TEST use random permutation
    for i=1:numel(Y)
        Y{i} = Y{i}(:,randperm(n));
    end
end

MarkerSize = 3;
nb = figure; 
plot(Y{1}(1,:), Y{1}(2,:), 'ok','MarkerFaceColor', 'k', 'MarkerSize',MarkerSize),
hold on
plot(Y{2}(1,:), Y{2}(2,:), 'ok','MarkerFaceColor', 'k', 'MarkerSize',MarkerSize)
axis equal  
  
%% 1) Compute SW2 projection from 2 point-clouds
  
% parameters setting
d0 = 10;
X0 = Y{1}; % randn(2,n)*.1+repmat([5;0], [1,n]);
w = [1 1]; w = w/sum(w); % weights for projection
options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'grad'; % 'grad' or 'stochastic'

options.niter = 1e3;
options.eps = 1e-6;
options.display = 0;

SW2_bary = [];
SW2_bary_cost = [];
ntests = 3;
MarkerSize = 3;
for i=1:ntests
    
    if INCREASE_DIR && i>1
        d = d*10; 
    else
        d = d0; % 1e1; % number of directions
    end
    options.ndir = d;
    if RANDOM_DIR
        t = rand(1,d)*pi; 
        %t = t + rand*pi; % uniform sampling or orientations
    else
        t = (0:d-1)/d*pi; % uniform sampling or orientations
    end
    options.base = [cos(t); sin(t)]; % directions set
    
    if RANDOM_INIT
        X0 = randn(2,n)*1+3;
    end
    
    
    tic, [SW2_bary{i},E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w,options); toc
    switch mod(i,6)
        case 1, symbol_flag = 'o'; color_flag = 'r';
        case 2, symbol_flag = 'o'; color_flag = 'b';
        case 3, symbol_flag = 'o'; color_flag = 'g';
        case 4, symbol_flag = 'o'; color_flag = 'm';
        case 5, symbol_flag = 'o'; color_flag = 'y';
        otherwise, symbol_flag = 'o'; color_flag = 'c';
    end
    
    figure(nb),
    plot(SW2_bary{i}(1,:) , SW2_bary{i}(2,:) , symbol_flag, ... %'Color', color_flag, ...
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
    
    drawnow
    
    SW2_bary_cost(i) = 2*n*E(end)/d;
end

      
%% 2) comparison with exact OT (bipartite graph matching with Hungarian Algorithm or Mosek library)

if (n>Nmax_comp)
    disp('STOP : computation may be very long with such number of points !')
    
    return
else
    % build cost matrix : Cij = || Y1(i) - Y2(j) ||^2
    Y12 = repmat(sum(Y{1}.^2,1) , [n 1]); % L2^2 norm of Y1
    Y22 = repmat(sum(Y{2}.^2,1)', [1 n]); 
    costMat = Y12 + Y22 - 2*Y{2}'*Y{1}; % L2^2 distance : figure, imagesc(costMat)
    
    if (n<=1e2) % Hungarian Algorithm
        disp('Hungarian algorithm!')
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

W2_proj = Y{2}(:,perm);
W2_bary = 1/2*Y{1} + 1/2* W2_proj;
W2_bary_cost = 1/2*sum( sum(( W2_bary(:,:) - Y{1}(:,:) ).^2) + sum( ( Y{2}(:,perm)-W2_bary(:,:) ).^2 ) ) % display exact L2-OT cost

symbol_flag = 'o'; color_flag = 'k';
MarkerSize = 4;
figure(nb),
% plot(W2_bary(1,:), W2_bary(2,:), 'om')
plot(W2_bary(1,:), W2_bary(2,:), symbol_flag, ... 
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
%title('Wasserstein projection and its sliced approximation')
%legend('Y1', 'Y2', 'SW2_bary', 'W2_bary')
drawnow
axis off

figure, 
plot(SW2_bary_cost)


%% build figures for submission paper
MarkerSize = 3;
symbol_flag = 'o';

figure,
for i=1:2
    switch i
        case 1, X = Y{1}; color_flag = 'k';
        case 2, X = Y{2}; color_flag = 'k';
        case 3, X = SW2_bary{1}; color_flag = [0 0.4 1];
        case 4, X = SW2_bary{2}; color_flag = [1 0.4 0];
        case 5, X = W2_bary; color_flag = 'k';
    end

    plot(X(1,:), X(2,:), symbol_flag, ... 
            'MarkerEdgeColor', color_flag, ...
            'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
    axis equal, axis off
    
    if i<2, pause, end 
end

figure, hold on
N = numel(SW2_bary)
c = jet(N);
for i=1:N+1
    if i<=N
        X = SW2_bary{i}; color_flag = c(i,:);
    else
        X = W2_bary; color_flag = 'k';
    end

    plot(X(1,:), X(2,:), symbol_flag, ... 
            'MarkerEdgeColor', color_flag, ...
            'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
    axis equal, axis off
end


