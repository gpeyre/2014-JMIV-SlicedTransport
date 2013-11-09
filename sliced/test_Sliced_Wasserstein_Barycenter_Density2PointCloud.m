%% script to test the sliced wasserstein barycenter for point clouds sampling of pictures
% julien.rabin@unicaen.fr - 2013

clear all, close all, clc

save_figure = 1

%% 0) Load pictures

N = 2^9; % N^2 is the number of pixels

rep = 'images/';
filename = {
    'autruche2.jpg' , ...   % 1
    'duck.png'      , ...   % 2
    'frog.png'      , ...   % 3
    'octopus.png'   , ...   % 4
    'octopus2.png'  , ...   % 5
    'rabbit.jpg'    , ...   % 6
    'seahorse.png'  , ...   % 7
    'square.png'    , ...   % 8
    'circle.png'    , ...   % 9
    'triangle.png'  , ...   % 10
    'ring.png'      , ...   % 11
};

% selec = [1 2 6] % ostrich - duck - rabbit
% selec = [8 9 10] % square - circle - triangle
selec = [2 9 11] % duck - circle - ring

for j=1:numel(selec)
    
    % d_Y{j} = imread([rep, filename{selec(j)}]);
    im{j} = load_image([rep, filename{selec(j)}(1:end-4)], N);
    
    % rescale in [0,1] and binarize
    im{j} = rescale(double(im{j}));
    im{j} = im{j}>0.5;
    
    if size(im{j},3) == 3
        im{j} = sum(im{j},3) / 3;
    end
end
% display
figure, 
for j=1:numel(im)
    subplot(1,numel(im),j), imshow(im{j}), axis equal, axis off
end
drawnow


%% 0-bis) Generate point-clouds

disp('sampling ...')

n = 1e4; % number of points
tic
for j=1:numel(im)
    r = (1:size(im{j},1));
    c = (1:size(im{j},2));
   
    res = 2;
    if 0 % slow code with row and column pdf sampling
        for i=1:n
            [x(i),y(i)] = pinky(c,r,im{j},res);
        end
        Y{j} = [x(:)';y(:)']/N;
        
    else % fast sampling using 1D pdf
        I = imresize(im{j}, res);
        Y{j} = point_sampling(I,n) / N;
    end
end
toc


% shift to obtain a triangle
l = 10;
Y{2}(1,:) = Y{2}(1,:) + l;
Y{3}(1,:) = Y{3}(1,:) + l/2;
Y{3}(2,:) = Y{3}(2,:) + l*sqrt(3)/2;


% display
nb = figure; 
plot(Y{1}(1,:), Y{1}(2,:), 'r.'), hold on
plot(Y{2}(1,:), Y{2}(2,:), 'g.'),
plot(Y{3}(1,:), Y{3}(2,:), 'b.'), axis ij, axis off
Legende = {'Y1'; 'Y2'; 'Y3'; };
drawnow

%% 1) Compute barycenter from 3 point-clouds

w = { [1 1 1]; ... % iso-barycenter
      [4 1 1]; [1 4 1]; [1 1 4]; ...
      [1 1 0]; ...
      [1 0 1]; ...
      [0 1 1]; ...
      %[-1 1 2]; [1 -1 2]; ... % extra-polation
    }; % weights list

d = 1e2; % number of directions
t = (1:d)/d*pi; % uniform sampling or orientations
options.ndir = d;
options.base = [cos(t); sin(t)]; % directions set

options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'grad'; % 'grad' or 'stochastic'
options.nsubdir = max(3,round(d/10));
  
options.display = 0;

  
tic, 
for j=1:numel(w)
    % set initial position
    [tmp, maxi] = max(w{j}); % return first maximum 
    X0 = Y{maxi};
    
    % normalize weights
    % w{j} = w{j}/sum(w{j}(:));
    
    [SW2_Bary{j},Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,Y,w{j},options);

    switch j
        case 1, symbol_flag = 'o'; color_flag = 'k';
        case 2, symbol_flag = 'x'; color_flag = 'c';
        case 3, symbol_flag = '+'; color_flag = 'm';
        case 4, symbol_flag = 'd'; color_flag = 'y';
        otherwise,  symbol_flag = '.'; color_flag = 'k';
    end
    color_flag = w{j}/sum(w{j}(:));
    
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

% convolution with gaussian kernel
N = 5e2;
s = N/sqrt(n);

nb = figure;
bbox_Y = [];
for j=1:numel(Y)
    [d_Y{j}, bbox_Y{j}] = density_estimate2(Y{j},N,s);
    c_Y{j} = mean(Y{j},2);
    
    im_Y{j} = 255 - uint8(d_Y{j}/max(d_Y{j}(:))  * 255);
    figure(nb), imagesc(im_Y{j}),  colormap gray, axis image, axis off,
    drawnow
end

bbox_Bary = [];
for j=1:numel(w)
    [d_SW2_Bary{j}, bbox_Bary{j}] = density_estimate2(SW2_Bary{j},N,s);
    c_Bary{j} = mean(SW2_Bary{j},2);
    
    im_SW2_Bary{j} = 255 - uint8(d_SW2_Bary{j}/max(d_SW2_Bary{j}(:))  * 255);
    figure(nb), imagesc(im_SW2_Bary{j}),  colormap gray, axis image, axis off,
    drawnow
end

disp('ToDo : normalize intensity')

%% 3) build horizontal picture

% im = im_Y{1};
% for j = 2:numel(Y)
%     im = cat(2,im, im_Y{j} );
% end
% for j = 1:numel(w)
%     im = cat(2,im,im_SW2_Bary{j});
% end
% 
% if save_figure
%     imwrite(im,'results/barycenters_images.png');
%     ! open 'results/barycenters_images.png'
% else
%     figure, imagesc(im), colormap gray, axis off, axis equal
% end

%% 4) build big gray density map of all pictures
% disp('Warning : might be very slow ! Continue ?'), pause

M = 1e3;
N = round(M * sqrt(3)/2); % ratio length / width
if 0
    Ans = questdlg('Warning : might be very slow ! Continue ?', ...
                    'Picture Size :', ...
                    'No', '1e3', '5e3', '1e3');
    if strcmpi(Ans,'No')
        return
    else
       N = str2num(Ans) 
    end
end

s = sqrt(N*M)/sqrt(n)*1/3/2; % rule of thumb for uniform distribution
YY = cat(2, Y{:}, SW2_Bary{:});
[d_Bary, bbox] = density_estimate2(YY,[N M],s);


if save_figure
    im = 255 - uint8(d_Bary/max(d_Bary(:))  * 255); 
    imwrite(im,'results/barycenters_images.png');
    ! open 'results/barycenters_images.png'
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

s = sqrt(N*M)/sqrt(n)*1/3/3; % rule of thumb for uniform distribution
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
    imwrite(im,'results/barycenters_images_color.png');
    ! open 'results/barycenters_images_color.png'
else
    figure(nb), imshow(im),  colormap gray, axis image, axis off,
end
