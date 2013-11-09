% Geometric small figure to illustrate Sliced Wassterstein barycenters
% local minimum

close all
clear all
clc

N = 1; % Number of points
N = N*3*4*5; % required to get polygone figures
X = {};
fig_num = figure();
MarkerSize = 5;


a = 108/180*pi; % inner pentagone angle
b = pi-a; cb = cos(b); sb = sin(b);
c = b - (pi/2 - b); cc = cos(c); sc = sin(c);
L = 5;

% circle
t = (1:N)/N*2*pi;
X{1} = [cos(t); sin(t)]; % in [-1,1]^2
figure(fig_num), hold on,
plot(X{1}(1,:), X{1}(2,:), 'ro-', 'MarkerSize',MarkerSize)

% square
n = N/4;
t = (0:n-1)/n;
z = zeros(1,n);
o = z+1;
X{2} = [ [t;z] , [o;t] , [1-t;o] , [z;1-t] ]*2-1; % in [-1,1]^2
X{2}(1,:) = X{2}(1,:) + L;
figure(fig_num), hold on,
plot(X{2}(1,:), X{2}(2,:), 'mo-', 'MarkerSize',MarkerSize)

% equilateral triangle
n = N/3;
t = (0:n-1)/n;
z = zeros(1,n);

X{3} = [ [t;z] , [(2-t)/2 ; t*sqrt(3)/2] , [(1-t)/2 ; (1-t)*sqrt(3)/2] ];
X{3} = X{3}*2+repmat([-1;-sqrt(3)/2],[1,N]); % in [-1,1]^2
X{3}(1,:) = X{3}(1,:) + L*(1+cb);
X{3}(2,:) = X{3}(2,:) + L*sb;
figure(fig_num), hold on,
plot(X{3}(1,:), X{3}(2,:), 'bo-', 'MarkerSize',MarkerSize)

% pentagone
n = N/5;
t = (0:n-1)/n;
z = zeros(1,n);

X{4} = [ [(t-1)*cb;(1-t)*sb] , [t;z] , [1+t*cb ; t*sb] , [1+cb - t*sc ; sb + cc*t] , [1+cb-sc - t*sc ; sb + cc *(1-t)] ];
X{4} = X{4} + repmat([-1/2;-sb],[1,N]); % in [-1,1]^2
X{4}(1,:) = X{4}(1,:) + L*(1+cb-sc);
X{4}(2,:) = X{4}(2,:) + L*(sb+cc);
figure(fig_num), hold on,
plot(X{4}(1,:), X{4}(2,:), 'go-', 'MarkerSize',MarkerSize)

% star
o = z+1;
X{5} = [ [t*(1+cb);t*sb] , [1+cb-t*(1+2*cb);o*sb] , [t*(1+cb)-cb;(1-t)*sb] , [1 + t*(cb-sc) ; (sb + cc)*t] , [(1+cb-sc) * (1 - t); (sb + cc) *(1-t)] ];
X{5} = X{5} + repmat([-1/2;-sb],[1,N]); % in [-1,1]^2
X{5}(1,:) = X{5}(1,:) + L*(-cb);
X{5}(2,:) = X{5}(2,:) + L*sb;
figure(fig_num), hold on,
plot(X{5}(1,:), X{5}(2,:), 'yo-', 'MarkerSize',MarkerSize)

drawnow

%% compute barycenters from différent starting point

% paramaters
d = 1e3; % number of directions
t = (1:d)/d*pi + rand; % uniform sampling or orientations with random offset
options.ndir = d;

options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'grad';
options.niter = 1e3;
options.eps = 1e-14;
  
options.display = 0;

options.base = [cos(t); sin(t)]; % directions set

tic
for j=1:5
    [SW2_Bary{j},Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X{j},X,[1 1 1],options);
    switch j
        case 1, symbol_flag = 'o'; color_flag = 'r';
        case 2, symbol_flag = 'o'; color_flag = 'm';
        case 3, symbol_flag = 'o'; color_flag = 'b';
        case 4, symbol_flag = 'o'; color_flag = 'g';
        otherwise,  symbol_flag = 'o'; color_flag = 'y';
    end
    
    % display
    figure(fig_num), 
    plot(SW2_Bary{j}(1,:) , SW2_Bary{j}(2,:) , symbol_flag, ... %'Color', color_flag, ...
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
    drawnow
end
toc

return

%% figure with transparency
MarkerSize = 10;

t= 0:pi/10:2*pi;
r = .1;

figure, hold on
for j=1:5
    switch j
        case 1, symbol_flag = 'o'; color_flag = 'r';
        case 2, symbol_flag = 'o'; color_flag = 'm';
        case 3, symbol_flag = 'o'; color_flag = 'b';
        case 4, symbol_flag = 'o'; color_flag = 'g';
        otherwise,  symbol_flag = 'o'; color_flag = 'y';
    end
    
    plot(X{j}(1,:) , X{j}(2,:) , symbol_flag, ... %'Color', color_flag, ...
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',MarkerSize);
    
    for i=1:N
        pb=patch(SW2_Bary{j}(1,i) + r*cos(t),SW2_Bary{j}(2,i) + r*sin(t),color_flag,'edgecolor','none'); %,'facecolor',
        alpha(pb,.3);
    end
end
axis equal
axis off

