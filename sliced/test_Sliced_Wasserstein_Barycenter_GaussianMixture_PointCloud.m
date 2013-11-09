%% script to test the sliced wasserstein barycenter for point clouds sampling of gaussians mixture
% julien.rabin@unicaen.fr - 2013

clear all, close all, clc

save_figure = 1

%% 0) Generate gaussian mixture point-clouds

n = 1e3; % number of points PER GAUSSIAN
n_mix = 2; % number of gaussians in the mixture
N = n*n_mix; % total number of points
K = 3; % number of point clouds

% random parameters
seed = 1; % 1
rng(seed,'v5uniform')
rng(seed,'v5normal')

%theta = rand(K,n_mix)*pi; % angle
theta = repmat([0 1 2]'*(2*pi/3) + pi/6 + pi/2, [1,2]);
mu_x = repmat([-1,1], [3,1]) .* cos(theta); %3*(rand(K,n_mix)-.5); % mean
mu_y = repmat([-1,1], [3,1]) .* sin(theta); % 3*(rand(K,n_mix)-.5);
sx = .2*ones(K,n_mix) % 1/4 + 1/4*rand(K,n_mix); % first standard deviation
sy = .2*ones(K,n_mix) % 1/8 + 1/8*rand(K,n_mix); % second standard deviation

% generate point clouds
Y = []; Y{K} = [];
for j = 1:K
    for i = 1:n_mix
        x = randn(1,n) * sx(j,i);
        y = randn(1,n) * sy(j,i);
        t = theta(j,i);
        c = cos(t); s = sin(t);
        Y{j} = cat(2, Y{j}, [ mu_x(j,i) + c*x-s*y; mu_y(j,i) + s*x+c*y ]);
    end
end


% shift to obtain a triangle
l = 10;
Y{2}(1,:) = Y{2}(1,:) + l;
Y{3}(1,:) = Y{3}(1,:) + l/2;
Y{3}(2,:) = Y{3}(2,:) + l*sqrt(3)/2;


% display
nb = figure; 
plot(Y{1}(1,:), Y{1}(2,:), 'r.'), hold on
plot(Y{2}(1,:), Y{2}(2,:), 'g.'),
plot(Y{3}(1,:), Y{3}(2,:), 'b.'), axis off
Legende = {'Y1'; 'Y2'; 'Y3'; };
axis equal
% grid on, axis on
drawnow

%% 1) Compute barycenter from 3 point-clouds

w = { [1 1 1]; ... % iso-barycenter
      [3 1 1]; [1 3 1]; [1 1 3]; ...
      [3 1 0]; [1 3 0]; [0 1 3]; [3 0 1]; [0 3 1]; [1 0 3]; ...
      [1 1 0]; ...
      [1 0 1]; ...
      [0 1 1]; ...
    }; % weights list


d = 1e2; % number of directions
t = (1:d)/d*pi; % uniform sampling or orientations
options.ndir = d;

options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'grad'; % 'grad' or 'stochastic'
options.nsubdir = max(3,round(d/10));
options.eps = 1e-4; % stopping criterion

options.base = [cos(t); sin(t)]; % directions set
  
tic, 
for j=1:numel(w)
    % normalize weights
    p = w{j}/sum(w{j}(:));
    
    % set initial position
    if 0 % NN initialisation
        [tmp, maxi] = max(p); % return first maximum 
        X0 = Y{maxi};
        
    elseif 0 % random : isotropic gaussian pdf in the center of the weighted barycenter of the centers
        X0 = randn(2,N) + repmat( mean(cat(2,p(1)*Y{1},p(2)*Y{2},p(3)*Y{3}),2) , [1 N]);
        
    else % weighted average and subsampling of inputs 
        n1 = round(N*p(1));
        n2 = round(N*p(2));
        n3 = (N-n1-n2);
        ech = randperm(N);
        X0 =  cat(2, Y{1}(:,ech(1:n1)), Y{2}(:,ech(1:n2)), Y{3}(:,ech(1:n3)) );
        X0 = X0 + repmat( mean(cat(2,p(1)*Y{1},p(2)*Y{2},p(3)*Y{3}),2) - mean(X0,2) , [1 N]); % center
    end
    
    [SW2_Bary{j},Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w{j},options);

    switch j
        case 1, symbol_flag = 'o'; color_flag = 'k';
        case 2, symbol_flag = 'x'; color_flag = 'c';
        case 3, symbol_flag = '+'; color_flag = 'm';
        case 4, symbol_flag = 'd'; color_flag = 'y';
        otherwise,  symbol_flag = '.'; color_flag = 'k';
    end
    color_flag = p;
    
    % display
    figure(nb), 
    plot(SW2_Bary{j}(1,:) , SW2_Bary{j}(2,:) , symbol_flag, ...
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',2);
    
    Legende{end+1} = ['weigths : ' num2str(w{j})];
    drawnow
end
toc
legend(Legende)
      
      
%% 2) estimate 2D density

% convolution with isotropic gaussian kernel
N = 2^8; % N^2 is the number of pixels
s = N/100*4; % scale parameter of parzen window

nb = figure;
for j = 1:numel(Y)
    d_Y{j} = 255 - uint8( density_estimate(Y{j},N,s) * 255); % density
    c_Y{j} = mean(Y{j},2);
    
    figure(nb), imagesc(d_Y{j}), colormap gray, axis image, axis off,
    drawnow
end
for j=1:numel(w)
    d_SW2_Bary{j} = 255 - uint8( density_estimate(SW2_Bary{j},N,s) * 255);
    c_Y{j} = mean(SW2_Bary{j},2);
    
    figure(nb), imagesc(d_SW2_Bary{j}),  colormap gray, axis image, axis off,
    drawnow
end

disp('ToDo : normalize intensity')

%% 3) build picture

% im = d_Y{1} ;
% for j = 2:numel(Y)
%     im = cat(2,im, d_Y{j} );
% end
% for j = 1:numel(w)
%     im = cat(2,im, d_SW2_Bary{j});
% end
% 
% if save_figure
%     imwrite(im,'results/barycenters_Gaussian_Mixture.png');
%     ! open 'results/barycenters_Gaussian_Mixture.png'
% else
%     figure(nb), imagesc(im), colormap gray, axis equal, axis off
% end

%% 4) build big map of all pictures
% disp('Warning : might be very slow ! Continue ?'), pause

M = 1e3;
if 1
    Ans = questdlg('Warning : might be very slow ! Continue ?', ...
                    'Picture Size :', ...
                    'No', '1e3', '5e3', '1e3');
    if strcmpi(Ans,'No')
        return
    else
       M = str2num(Ans) 
    end
end
N = round(M * sqrt(3)/2); % ratio length / width

s = N/sqrt(n)*.8;
YY = cat(2, Y{:}, SW2_Bary{:});
[d_Bary, bbox] = density_estimate2(YY,[N M],s);

if save_figure
    im = 255 - uint8(d_Bary/max(d_Bary(:))  * 255); 
    imwrite(im,'results/barycenters_gaussians_Mixture.png');
    ! open 'results/barycenters_gaussians_Mixture.png'
else
    figure(nb), imagesc(d_Bary),  colormap gray, axis image, axis off,
end

%% 5) build big color density map of all pictures
% disp('Warning : might be very slow ! Continue ?'), pause

M = 1e3;
if 1
    Ans = questdlg('Warning : might be very slow ! Continue ?', ...
                    'Picture Size :', ...
                    'No', '1e3', '5e3', '1e3');
    if strcmpi(Ans,'No')
        return
    else
       M = str2num(Ans) 
    end
end
N = round(M * sqrt(3)/2); % ratio length / width

s = sqrt(N*M)/sqrt(n)*1/3; % rule of thumb for uniform distribution
Z = zeros(N,M,1);
d_Bary = zeros(N,M,3);
% bbox = [0,l*sqrt(3)/2;0,l];
for i=1:numel(Y)
    p = zeros(3,1); p(i) = 1;
    tmp = density_estimate2(Y{i},[N M],s,bbox);
    d_Bary = d_Bary + cat(3,p(1)*tmp,p(2)*tmp,p(3)*tmp);
end
for i=1:numel(SW2_Bary)
    p = w{i}/sum(w{i}(:));
    tmp = density_estimate2(SW2_Bary{i},[N M],s,bbox);
    d_Bary = d_Bary + cat(3,p(1)*tmp,p(2)*tmp,p(3)*tmp);
end

im = 255 - uint8(d_Bary/max(d_Bary(:))  * 255); 
if save_figure
    imwrite(im,'results/barycenters_gaussians_Mixture_color.png');
    ! open 'results/barycenters_gaussians_Mixture_color.png'
else
    figure(nb), imshow(im),  colormap gray, axis image, axis off,
end
