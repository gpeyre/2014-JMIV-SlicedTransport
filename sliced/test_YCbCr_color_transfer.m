
clear all, close all, clc  

%% Load 2 images

if 0 % rainbow color transfer
    im{1} = imread('peppers.png');
    
    % creer arc-en-ciel
    row = 50;
    col = 300;
    H = repmat( linspace(0,1,col), [row, 1]);
    rgb = hsv2rgb( cat(3,H, ones(row,col,2)) );
    im{2} = uint8(rgb*255);
end


N = 2^8; % N^2 is the number of pixels
N2 = N*N;

rep = 'images/';
filename = {
    'colors-1.jpg' , ...   % 1
    'colors-2.jpg' , ...   % 2
    'fleur-1.jpg'  , ...   % 3
    'fleur-2.jpg'  , ...   % 4
    'parrot-1.jpg' , ...   % 5
    'parrot-2.jpg' , ...   % 6
    'wheat-1.jpg'  , ...   % 7
    'wheat-2.jpg'  , ...   % 8
};

selec = [5 6];
K = numel(selec);

im = [];
for i=1:K
    % load in rgb
    im{i} = uint8(load_image([rep, filename{selec(i)}(1:end-4)], N));
    
    % conversion in YCbCr
    YCbCr{i} = rgb2ycbcr(im{i});
    
    % separate luminance channel 'Y' from chroma
    Y{i}  = YCbCr{i}(:,:,1); % in [16 235]/255
    Cb{i} = YCbCr{i}(:,:,2); % in [16 240]/255
    Cr{i} = YCbCr{i}(:,:,3); % in [16 240]/255
    CbCr{i} = [Cb{i}(:), Cr{i}(:)]'; % [2 x N] point cloud
    
    % compute chroma histograms
    q_Cr = (16:1:240);
    CbCr_hist{i} = hist3(CbCr{i}',{16:1:240 16:1:240})/N/N;
end

% display pictures and chroma
figure, 
for i=1:K
    subplot(3,K,i), imshow(im{i}), axis equal, axis off, title(['im' num2str(i)])
    Chroma{i} = ycbcr2rgb( cat(3, 127*ones(N), Cb{i}, Cr{i}) ); %  set for luminance Y=127/255
    subplot(3,K,i+K), imshow(Chroma{i}), axis equal, axis off, title(['Chroma ' num2str(i)])
    subplot(3,K,i+2*K), imagesc(CbCr_hist{i}), axis equal, axis off, title(['Pdf Chroma ' num2str(i)])
end
drawnow

%% Visualisation nuage de point
 
% nuage de point
figure,
plot(Cb{1}(:), Cr{1}(:), 'b.'), hold on,
plot(Cb{2}(:), Cr{2}(:), 'm.'), 
xlabel('Cb'), ylabel('Cr')

% nuage de point en couleur
for i=1:K
    Color{i} = cat(3,ones(N)*127,Cb{i},Cr{i});
    
    Color{i} = double(ycbcr2rgb(uint8(Color{i}))); % figure, imshow(Color)
    Color{i} = reshape(Color{i}, [N2 3]);
    Nech = 1e3;
    ech = randperm(N2); ech = ech(1:Nech);
    figure,
    scatter(Cb{i}(ech),Cr{i}(ech), 50,  Color{i}(ech,:)/255, 'fill') %
    xlabel('Cb'), ylabel('Cr')
end

%% transfert de couleur avec Sliced Wasserstein energy on point clouds

% parameters setting
X0 = double(CbCr{1});
X1 = double(CbCr{2});
w = [0 1]; w = w/sum(w); % weights for projection
d = 10; % number of directions
t = (1:d)/d*pi; % uniform sampling or orientations
options.ndir = d;
options.base = [cos(t); sin(t)]; % directions set
options.step = 1; % descent step (maximum)
options.hessian = 1; % use hessian normalization
options.method = 'grad'; % 'grad' or 'stochastic' or 'average'
options.nsubdir = round(d/10);

options.display = 0;

tic, [SW2_proj,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,{X0,X1},w,options); toc

% image transférée
im12 = ycbcr2rgb( cat(3,Y{1},reshape(SW2_proj',[N N 2])) );
figure, 
imshow(im12)

%% transfert de couleur avec Sliced Wasserstein energy on density (weighted dirac distributions)

error('A FINIR')

addpath('../radon');
    
use_fast_slant_stack = 0; % matlab's radon or not

N_imgs = 2; % number of pdfs

% precompute Radon transforms
tic;
radon_pdfs = cell(N_imgs,1);
if (use_fast_slant_stack)
    addpath('../fast-slant-stack/SlantStackP');
    for i=1:N_imgs
        radon_pdfs{i} = abs(real(FastSlantStack(CbCr_hist{i})))*N*2/sqrt(2);
    end;
else
    for i=1:N_imgs
        radon_pdfs{i} = radon(CbCr_hist{i}, 0:179); % changer ces valeurs pour passer à 255 ou à N ?
    end
end
toc;

% compute radon histogram matching
N_interp = 1;
tic;
result = zeros( N, N, N_interp);
% matlabpool(8);
parfor i=1:N_interp
    
    radonInterp = interp_pdfs_Nd_radon(radon_pdfs{:}, 0, 1);
    
    if (use_fast_slant_stack)
        result(:,:,i) = real(Inv_FastSlantStack(radonInterp))*sqrt(2)/(N*2);
    else
        result(:,:,i) = iradon(radonInterp, 0:179, 'linear', 'Hann', 1., N);
    end;
    
end;
matlabpool close;
toc

return
%% Image du plan CbCr pour différentes intensités lumineuses
y = {16,127,235} % in [16 235]

x = 16:240;
[X,Y] = meshgrid(x);
figure
for i=1:numel(y)
   CbCr_RGB{i} = ycbcr2rgb( uint8(cat(3,y{i}*ones(size(X)), X, Y)) ); 
   subplot(1,3,i), imshow(uint8(CbCr_RGB{i}))
   title(['Y = ' num2str(y{i})])
end


