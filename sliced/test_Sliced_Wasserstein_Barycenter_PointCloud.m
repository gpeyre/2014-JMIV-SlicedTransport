%% script to test the sliced wasserstein barycenter for synthetic point clouds
% julien.rabin@unicaen.fr - 2013

clear all,
close all, 
%clc

save_figure = 1

%% 0) Generate discrete distributions
n = 2e3; % data size
K = 3; % number of input point-clouds

if 0 % generate periodic curves
    
    x = linspace(0,1,n)*2*pi;
    Y = { [cos(x); .5*sin(x) ]; ...
        [.5*cos(x); .8*sin(2*x) ]; ...
        [x/6.*cos(2*x); x/6.*sin(2*x) ]; ...
      };

    % Y{end+1} = [cos(x)+.1*cos(15*x); sin(x)+.1*sin(15*x) ];  % flower
    % Y{end+1} = [cos(x); .5*sin(x); ]; % ellipse
    % Y{end+1} = [.5*cos(x)+2; sin(x); ]; % ellipse
    % Y{end+1} = [ 1.2*sin(2*x); .5*cos(x); ]; % 8
    % Y{end+1} = [.5*cos(x); 1.2*sin(2*x); ]; % infinity
    % Y{end+1} = [cos(2*x).^3+4; sin(3*x).^3; ]; % arrow
  
elseif 1 % segments
    
    theta = [1/6; 1/6+1/2; 1/6+1/10]*pi; % angle
    si = [3, 3, 6]; % 1/4 + 1/4*rand(K,n_mix); % first standard deviation
   
    % generate point clouds
    Y = []; Y{3} = [];
    for j = 1:3
            %x = linspace(-1/2,1/2,n);
            x = rand(1,n)-.5;
            
            x = x * si(j);
            y = (rand(1,n)-.5)*.1;
            t = theta(j);
            c = cos(t); s = sin(t);
            Y{j} = [ c*x-s*y; s*x+c*y ];
    end

else % mixture of two gaussians
    n1 = round(n/4);
    n2 = n-n1;

    theta = [1/3,0; 2/3,1/4; 1/2,1]*pi % (0:K-1)/K*pi + pi/3; % angle
    sx = [1,1; 1,2; 2,1]/4; % first standard deviation
    sy = [1,1; 2,1; 1,1]/4; % second standard deviation
    mux = [1,-1; -1,1; -1,1];
    muy = [1,-1; 2,-1; 1,-1];
    
    for j = 1:K
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
l = 20; % distance between exemplar distributions centers
Shift = [0,0; l,0; l/2,l*sqrt(3)/2]'; % [2 x K] shifting parameters
use_shift = 0; % if 1, translate data for computing, else only for figure
if use_shift % translate input values
    for i=1:K
        Y{i}(1,:) = Y{i}(1,:) + Shift(1,i);
        Y{i}(2,:) = Y{i}(2,:) + Shift(2,i);
    end
    Shift_Y = zeros(2,K);
else
    Shift_Y = Shift;
end

% display
nb = figure; 
plot(Y{1}(1,:) + Shift_Y(1,1), Y{1}(2,:) + Shift_Y(2,1), 'r.'), hold on
plot(Y{2}(1,:) + Shift_Y(1,2), Y{2}(2,:) + Shift_Y(2,2), 'g.'),
plot(Y{3}(1,:) + Shift_Y(1,3), Y{3}(2,:) + Shift_Y(2,3), 'b.'), axis off
Legende = {'Y1'; 'Y2'; 'Y3'; };

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
options.nsubdir = round(d/10);
options.eps = 1e-4;
  
options.display = 0;

options.base = [cos(t); sin(t)]; % directions set

Shift_Bary = zeros(2,numel(w));


tic, 
for j=1:numel(w)
    % normalize weights
    p = w{j}/sum(w{j}(:));
    
    % set initial position
    if 0 % NN initialisation
        [tmp, maxi] = max(p); % return first maximum 
        X0 = Y{maxi};
        
    elseif 1 % random : isotropic gaussian pdf in the center of the weighted barycenter of the centers
        X0 = 1e-3*randn(2,n) + repmat( mean(cat(2,p(1)*Y{1},p(2)*Y{2},p(3)*Y{3}),2) , [1 n]);
        
    else % weighted average and subsampling of inputs 
        n1 = round(n*p(1));
        n2 = round(n*p(2));
        n3 = (n-n1-n2);
        ech = randperm(n);
        X0 =  cat(2, Y{1}(:,ech(1:n1)), Y{2}(:,ech(1:n2)), Y{3}(:,ech(1:n3)) );
        X0 = X0 + repmat( mean(cat(2,p(1)*Y{1},p(2)*Y{2},p(3)*Y{3}),2) - mean(X0,2) , [1 n]); % center
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
    
    if ~use_shift % need to shift barycenter for figure using shifting property of the barycenter
        Shift_Bary(:,j) = sum( repmat(p,[2 1]) .* Shift, 2);
    end
    
    % display
    figure(nb), 
    plot(SW2_Bary{j}(1,:) + Shift_Bary(1,j) , SW2_Bary{j}(2,:) + Shift_Bary(2,j) , symbol_flag, ...
        'MarkerEdgeColor', color_flag, ...
        'MarkerFaceColor', color_flag, 'MarkerSize',2);
    
    Legende{end+1} = ['weigths : ' num2str(w{j})];
    drawnow
    
    %pause
end
toc
legend(Legende)
      
      
%% 2) estimate 2D density

% convolution with gaussian kernel
N = 2^10; % N^2 : number of pixels
s = N/100*2;

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

% im = d_Y{1};
% for j = 2:numel(Y)
%     im = cat(2,im,d_Y{j});
% end
% for j = 1:numel(w)
%     im = cat(2,im,d_SW2_Bary{j});
% end
% 
% if save_figure
%     imwrite(im,'results/barycenters.png');
%     ! open 'results/barycenters.png'
% else
%     figure(nb), imagesc(im), colormap gray, axis off, axis equal
% end

%% 4) build big map of all pictures
% disp('Warning : might be very slow ! Continue ?'), pause

N = 1e3;
if 1
    Ans = questdlg('Warning : might be very slow ! Continue ?', ...
                    'Picture Size :', ...
                    'No', '1e3', '5e3', '1e3');
    if strcmpi(Ans,'No')
        return
    else
       N = str2num(Ans) 
    end
end

s = N/sqrt(n)*.1;
YY = [];
% shifting all point clouds for figure
for i=1:numel(Y)
    YY = cat(2, YY, Y{i} + repmat(Shift_Y(:,i),[1 n]) );
end
for i=1:numel(SW2_Bary)
    YY = cat(2, YY, SW2_Bary{i} + repmat(Shift_Bary(:,i),[1 n]) );
end

[d_Bary, bbox] = density_estimate2(YY,N,s);

if save_figure
    im = 255 - uint8(d_Bary/max(d_Bary(:))  * 255); 
    imwrite(im,'results/barycenters.png');
    ! open 'results/barycenters.png'
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
    tmp = density_estimate2(Y{i} + repmat(Shift_Y(:,i),[1 n]),[N M],s,bbox);
    d_Bary = d_Bary + cat(3,p(1)*tmp,p(2)*tmp,p(3)*tmp);
end
for i=1:numel(SW2_Bary)
    p = w{i}/sum(w{i}(:));
    tmp = density_estimate2(SW2_Bary{i} + repmat(Shift_Bary(:,i),[1 n]) ,[N M],s,bbox);
    d_Bary = d_Bary + cat(3,p(1)*tmp,p(2)*tmp,p(3)*tmp);
end

im = 255 - uint8(d_Bary/max(d_Bary(:))  * 255); 
if save_figure
    imwrite(im,'results/barycenters_color.png');
    ! open 'results/barycenters_color.png'
else
    figure(nb), imshow(im),  colormap gray, axis image, axis off,
end
