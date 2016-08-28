%% Test for point-cloud sampling using an input image as a pdf and Floyd regularization

close all, 
clear all,
clc

image_resize = 2^8; % if >0 number of rows

%% 1) load natural or synthetic image and convert to 2D pdf

if 1 % natural image
    
    filename = 'cameraman.tif' % 'cameraman.tif', 'board.tif', 'peppers.png', 'football.jpg'
    im = imread(filename);
    if size(im,3)>1
        disp('RGB -> grey level conversion')
        im = rgb2gray(im);
    end

    % picture resizing
    if image_resize
        im = imresize(im, [image_resize NaN]); % prescribe number of rows
    end
    
else % synthetic image
    
    n = 2^7;
    if 1 % ramp
        im = repmat((1:n), [n/8 1]);
        
    elseif 0 % gaussians
        gauss = @(x,mu,sig) exp( - (x-mu).^2 /(2*sig^2) );
        gauss2 = @(X,Y,mu,sig) gauss(X,mu(1),sig(1)) * gauss(Y,mu(2),sig(2));

        [Y,X] = meshgrid(1:n);
        im = 3*gauss2(X,Y,[.2 .3]*n,[1 .5]*n/10) + gauss2(X,Y,[.8 .7]*n,[1/3 1]*n/8);
    elseif 0 % square
        im = zeros(n);
        im(n/4:n/4*3, n/4:n/4*3) = 1;
        im(n*3/8:n/8*5, n*3/8:n/8*5) = 0;
        w = gausswin(10,4);
        im = conv2(im,w*w','same');
    end
end

[n,m,d] = size(im);
N = n*m;

im = double(im)/sum(im(:)); % picture processed as a pdf (sum of weighted dirac masses)

figure(1), imagesc(im), colormap gray, axis image, axis off

if 0 % TEST : display level-set using periodic contrast enhancement
    t = mod( (1:256)*4 , 256);% randperm(256)-1;
    figure(2), imshow(t(uint8(im*255/max(im(:)))+1), [0 255])
end

%% 2)  discretize pdf using random sampling

K = round( min(4e7/N, N/10) ) % number of centroids / points with unit mass
niter = 10;
res = 2; % up-sampling factor

% random centroids
r = (1:n); c = (1:m);
if 0 % with uniform distribution
    Ci = randi(n, 1, K); Ci = r(Ci);
    Cj = randi(m, 1, K); Cj = c(Cj);
    
elseif 1 % sampling using 1-D pdf of intensity
    I = imresize(im, res);
    C = point_sampling(I,K);
    Ci = C(2,:)/res; Cj = C(1,:)/res;
    
else % sampling using column and row pdf (nasty method : should use poisson disk algorithm instead for instance)
    for it=1:K
        [Cj(it),Ci(it)]=pinky(c,r,im,res);
    end
    
end
C = [Ci(:)'; Cj(:)']; % [2 x K] : figure, plot(C(2,:), C(1,:), '.')

fig_num = figure,
    imagesc(im), colormap gray, axis image, axis off, hold on, 
    scatter(C(2,:), C(1,:), 50, [0 0 .5], 'fill'); 
    scatter(C(2,:), C(1,:), 30, [1 1 0],  'fill'); hold off


%% 3) Regularize point cloud with Loyd algorithm

% precomputation
r = (1:n); c = (1:m);
[J,I] = meshgrid(c,r); % pixels coordinates : figure, imagesc(I)
X = [I(:)'; J(:)']; % [2 x N] : figure, plot(X(1,:), X(2,:), '.')
X2 = repmat(sum(X.^2,1),[K 1]); % square L2 norm matrix [K x N]
V = im(:)'; % intensity values (pdf of weighted dirac sum in [0,1])
V2 = cat(1,V,V);
VV = repmat(V,[K 1]); % [K x N] 
XX = repmat(reshape(X', [1 N 2]), [K 1]); % [K x N x 2]

C2 = repmat(sum(C.^2,1)',[1 N]); % [K x N]
Energ = [];

% algorithm
tic
for it = 1:niter
   % step 1 : update labelling by NN (L2 metric)
   D2 = X2 + C2 - 2 * C'*X; % [K x N] matrix of L_2^2 distances
   [Mini, Lab] = min(D2,[],1);
   
   % energy = average distorsion
   Mini = max(0, Mini);
   Energ(end+1) = sum(Mini);
   
   % step 2 : weighted barycenter (using pdf value)
   if 0 % code with inner loop
       for it2=1:K
           id = Lab==it2;
           num = sum(id);
           if num>0
               W(it2) = sum( V(id) , 2);
               C(:,it2) = sum( X(:,id) .* V2(:,id) , 2);
           else
               %new random position ? 
               W(it2) = 1;
               C(:,it2) = [randi(n); randi(m);];
           end
       end
   else % code without inner loop but with huge matrices
       id = repmat(Lab, [K 1]) == repmat((1:K)', [1 N]);
       W = VV .*id; % [K x N] weights for each cluster
       C = reshape(sum( XX .* cat(3,W,W) , 2), [K 2])'; % [K 1 2] -> [2 K]
       W = sum( W , 2)'; % [2 K]
       
       % traitement des clusters vides
       id = W==0;
       if sum(id)>0
           W(id) = 1;
           C(:,id) = randi(min(n,m),2,sum(id));
       end
   end
   C = C./ repmat( W, [2 1] );
   C2 = repmat(sum(C.^2,1)',[1 N]); % [K x N]
   
   
   % display
   if 1 && mod(it*10, niter) == 0
       figure(fig_num), imagesc(im), colormap gray, axis image, axis off
       hold on, scatter(C(2,:), C(1,:), 30, [1 1 0], 'fill'); hold off
       title(['iterate #' num2str(it)])
       drawnow
        %pause
   end
end
toc

% display energy
figure, plot(Energ), title('Distorsion')

% final display
figure(fig_num), 
    imagesc(im), colormap gray, axis image, axis off, hold on, 
    scatter(C(2,:), C(1,:), 50, [0 0 .5], 'fill'); 
    scatter(C(2,:), C(1,:), 30, [1 1 0],  'fill'); hold off
    
% figure, imagesc(reshape(Lab,[n m])) % clustering regions
% figure, imagesc(reshape(sqrt(max(0,Mini)),[n m])) % distance to nearest centroid
% figure, imagesc(reshape(W(Lab),[n m])) % total mass of the region (should be almost uniform)

u = im/max(im(:)); u = cat(3,u,u,u);
rgb = hsv(K);

% Color_Lab = arrayfun(@(x)(rgb(x)), Lab);
x = reshape(rgb(Lab,:),[n m 3]);
figure, imshow(x.*(1+2*u)/3)
hold on, 
scatter(C(2,:), C(1,:), 15, [0 0 .5], 'fill'); 
scatter(C(2,:), C(1,:), 10, [1 1 0],  'fill');
hold off

%% 2bis) Compute density from the point-cloud using parzen window
X = [C(2,:); C(1,:)];
im0 = density_estimate(X,size(im,1),5);

figure, imagesc(im0), colormap gray, axis image, axis off,

%% 3) alternative method : blue noise sampling or farthest point sampling

disp('todo')

