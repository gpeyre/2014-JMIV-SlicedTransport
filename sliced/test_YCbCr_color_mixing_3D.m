% Script to test color mixing & transfer using chromatic components AND
% luminance of images in YCbCr color space instead of the RGB coordinates
% This script uses also a post-processing regularization to avoid small
% artifacts

clear all, 
%close all, clc 
RNG_STATE = rng(0,'v5uniform');

im_save = 0
use_proj = 0 % permet de choisir si on veut projeter les images sur le barycentre ou recalculer le barycentre en initialisant de l'image 
perform_grading = 0

% build direction set
    d = 3; % number of values for latitude AND longitude
    disp(['Number of directions : ', num2str(d^2)])
    phi = (1:d)/d*2*pi; % uniform sampling of longitude
    psi = acos( (0:1:d-1)/d ); % uniform sampling or orientations
    [psi,phi] = meshgrid(psi,phi);
    phi = phi + repmat((0:d-1)/d*2*pi/d , [d,1] ); % offset in longitude
    psi = psi(:)';  phi = phi(:)';
    B = [cos(phi).*sin(psi); sin(phi).*sin(psi); cos(psi)];
    if 1, figure, plot3( B(1,:), B(2,:), B(3,:) ,'o'), end
    
    
% SW2 parameters
    options.ndir = d*d;
    options.base = B; % directions set
    options.step = 1; % descent step (maximum)
    options.hessian = 1; % use hessian normalization
    options.method = 'grad'; % 'grad' or 'stochastic' or 'average'
    %options.nsubdir = max(3,round(d/10));

    options.display = 0;


% ToDo
sprintf( ...
    ['ToDo list :\n', ...
    '- try other color space (Lab)\n', ...
    '- interpolate chromatic OT map with Nicolas code\n', ...
    '- application to eigen face for average template\n', ...
    ])

%% Load 2 images
disp(' '), disp('Load images ...'), disp(' ')
if 0 % rainbow color transfer
    im{1} = imread('peppers.png');
    
    % creer arc-en-ciel
    row = 50;
    col = 300;
    H = repmat( linspace(0,1,col), [row, 1]);
    rgb = hsv2rgb( cat(3,H, ones(row,col,2)) );
    im{2} = uint8(rgb*255);
end


rep = 'images/';
filename = {
    'colors-1.jpg' , ...   % 1
    'colors-2.jpg' , ...   % 2
    'fleur-1.jpg'  , ...   % 3
    'fleur-2.jpg'  , ...   % 4
    'fleur-3.jpg'  , ...   % 5
    'parrot-1.jpg' , ...   % 6
    'parrot-2.jpg' , ...   % 7
    'parrot-3.jpg' , ...   % 8
    'wheat-1.jpg'  , ...   % 9
    'wheat-2.jpg'  , ...   % 10
    'Notre_Dame_des_flots_1.jpg'  , ...   % 11
    'Notre_Dame_des_flots_2.jpg'  , ...   % 12
    'Notre_Dame_des_flots_3.jpg'  , ...   % 13
    'Van_Gogh_1.jpg'  , ...   % 14
    'Van_Gogh_2.jpg'  , ...   % 15
    'Van_Gogh_3.jpg'  , ...   % 16
    'Van_Gogh-4.jpg'  , ...   % 17
    'Van-Gogh_5.jpg'  , ...   % 18
    'Van-Gogh_6.jpg'  , ...   % 19
    'Van-Gogh-7.jpg'  , ...   % 20
    'clockHD-1.jpg'   , ...   % 21
    'clockHD-2.jpg'   , ...   % 22
    'clockHD-3.jpg'   , ...   % 23
    'clockmontague-1.jpg'  , ...   % 24
    'clockmontague-2.jpg'  , ...   % 25
    'clockmontague-3.jpg'  , ...   % 26
    '4625786629_1808574482_b.jpg'  , ...   % 27
    '8245632244_78ca92ed42_h.jpg'  , ...   % 28
    '8581254643_8ead330c4c_h.jpg'  , ...   % 29
};

selec = [24 25 26]; % 6 7 8 % 3 4 5 % 11 12 13 % 14 15 17 % 6 7 14 % 27 28 29
K = numel(selec);
if (K<2 || K>3)
    error('Error : number of pictures should be 2 ou 3 only')
end



N = 2^9; % N^2 is the number of pixels
im = load_image([rep, filename{selec(1)}(1:end-4)]); % ,N
[n m c] = size(im);
M = round(N*m/n);
if mod(M,2)==1,
    M=M+1;
end
N2 = N*M;

im = [];
for i=1:K
    % load in double
    im{i} = double(uint8(load_image([rep, filename{selec(i)}(1:end-4)]) )); %  % uint8 -> double cast values in [0,255]
    %[n m c] = size(im{i});
    im{i} = image_resize(im{i},N,M,c);
    %im{i} = crop(im{i}, N, [N N]/2)
    
    % conversion in YCbCr
    YCbCr{i} = rgb2ycbcr(im{i});
    
    % separate luminance channel 'Y' from chroma
    Y{i}  = YCbCr{i}(:,:,1); % in [16 235]
    Cb{i} = YCbCr{i}(:,:,2) + 128; % in [16 240] if adding +128 for double format (BUG of matlab 'rgb2ycbcr' function, unnecessary for uint8 format)
    Cr{i} = YCbCr{i}(:,:,3) + 128; % in [16 240]
    YCbCr{i} = [Y{i}(:), Cb{i}(:), Cr{i}(:)]'; % [3 x N] point-cloud
    
    % compute chroma histograms
    q_Cr = (16:1:240);
    CbCr_hist{i} = hist3(YCbCr{i}(2:3,:)',{16:1:240,16:1:240})/N/N;
end

% display pictures and chroma
figure(1), 
for i=1:K
    subplot(3,K,i), imshow(uint8(im{i})), axis equal, axis off, title(['im' num2str(i)])
    Chroma{i} = ycbcr2rgb( cat(3, 127*ones(N,M), Cb{i}, Cr{i})/ 255 ) * 255; %  set for luminance Y=127/255
    subplot(3,K,i+K), imshow(uint8(Chroma{i})), axis equal, axis off, title(['Chroma ' num2str(i)])
    subplot(3,K,i+2*K), imagesc(CbCr_hist{i}), axis equal, axis off, title(['Pdf Chroma ' num2str(i)])
end
drawnow

%% Visualisation en nuages de points
 
% nuage de point
figure(2),
plot3(Y{1}(:), Cb{1}(:), Cr{1}(:), 'b.'), hold on,
plot3(Y{2}(:), Cb{2}(:), Cr{2}(:), 'm.'), 
xlabel('Cb'), ylabel('Cr')
if K==3
    plot3(Y{3}(:), Cb{3}(:), Cr{3}(:), 'r.'), 
end
drawnow

% nuages de points en couleur
figure(1),
for i=1:K
    Color{i} = cat(3,ones(N,M)*127,Cb{i},Cr{i});
    
    Color{i} = double(ycbcr2rgb(uint8(Color{i}))); % figure, imshow(Color)
    Color{i} = reshape(Color{i}, [N2 3]);
    Nech = 1e3;
    ech = randperm(N2); ech = ech(1:Nech);
    
    subplot(3,K,i+2*K)
    scatter(Cb{i}(ech),Cr{i}(ech), 50,  Color{i}(ech,:)/255, 'fill') %
    xlabel('Cb'), ylabel('Cr')
end
drawnow

%% Step 1 : compute color average (and projection of image 1 if option 'use_proj')

disp(' '), disp('Compute Barycenter ...'), disp(' ')
if K==3
    X1 = double(YCbCr{1});
    X2 = double(YCbCr{2});
    X3 = double(YCbCr{3});
    w = [1 1 1]; w = w/sum(w); % weights for projection
   
    if use_proj % random init
        X0 = randn(size(X1));
    else % use image 1 as initialisation to compute barycenter AND ALSO the projection of image1 onto it
        X0 = X1;
    end
    tic, [SW2_bary0,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,{X1,X2,X3},w,options); toc

    if ~use_proj
        SW2_bary1 = SW2_bary0;
    
        % image transférée
        im1Bary = ycbcr2rgb( reshape(SW2_bary1',[N M 3]) / 255 ) * 255;
        figure(10), 
        subplot(2,K,1), imshow(im1Bary/255)
    end
    
    figure(2),
    plot3(SW2_bary0(1,:), SW2_bary0(2,:), SW2_bary0(3,:), 'gx'), 
    
else
    SW2_bary1 = double(YCbCr{2});
end
drawnow

%% Step 2 : transfert de couleur des autres images avec Sliced Wasserstein energy on point clouds

disp(' '), disp('Project images ...'), disp(' ')
if K==2
    % parameters setting
    X1 = double(YCbCr{1});
    X2 = double(YCbCr{2});
    X3 = double(YCbCr{3});
    
    w = [0 1]; w = w/sum(w); % weights for projection

    tic, [SW2_proj,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X0,{X0,X1},w,options); toc

    % image transférée
    im12 = ycbcr2rgb( reshape(SW2_proj',[N N 3]) / 255 ) * 255;
    figure(10), 
    imshow(im12/255)
else
    if ~use_proj % recompute barycenter from different initialisation
        
        tic, [SW2_bary2,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X2,{X1,X2,X3},w,options); toc
        
        tic, [SW2_bary3,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X3,{X1,X2,X3},w,options); toc
       
    else % if use_proj  use previous barycenter
        
        w = [0 1]; w = w/sum(w); % weights for projection
        
        options.method = 'stochastic'; % 'grad' or 'stochastic' 
        options.nsubdir = options.ndir; % max(3,round(d/10));
        
        tic, [SW2_bary1,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X1,{X1,SW2_bary0},w,options); toc
        
        tic, [SW2_bary2,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X2,{X1,SW2_bary0},w,options); toc
        
        tic, [SW2_bary3,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X3,{X1,SW2_bary0},w,options); toc
       
        
        im1Bary = ycbcr2rgb( reshape(SW2_bary1',[N M 3]) / 255 ) * 255;
        figure(10), subplot(2,K,1), imshow(im1Bary/255)
    end
    
    % nuages de points barycentre
    figure(2), plot3(SW2_bary2(1,:), SW2_bary2(2,:), SW2_bary2(3,:), 'y+'),
    
    % image transférée
    im2Bary = ycbcr2rgb( reshape(SW2_bary2',[N M 3]) / 255 ) * 255;
    figure(10), subplot(2,K,2), imshow(im2Bary/255)

    % nuage de point barycentre
    figure(2), plot3(SW2_bary3(1,:), SW2_bary3(2,:), SW2_bary3(3,:), 'co'), 
    % image transférée
    im3Bary = ycbcr2rgb( reshape(SW2_bary3',[N M 3]) / 255 ) * 255;
    figure(10), subplot(2,K,3), imshow(im3Bary/255)
    
    % nuages de points en couleur
    Color{1} = im1Bary; Color{2} = im2Bary; Color{3} = im3Bary;
    SW2_bary{1} = SW2_bary1; SW2_bary{2} = SW2_bary2; SW2_bary{3} = SW2_bary3;
    figure(10),
    for i=1:K
        Color{i} = reshape(Color{i}, [N2 3]);

        subplot(2,K,i+K)
        scatter(SW2_bary{i}(1,ech),SW2_bary{i}(2,ech), 50,  Color{i}(ech,:)/255, 'fill') %
        xlabel('Cb'), ylabel('Cr')
    end

end

if im_save && K==3
   rep = 'results/'
   
   if use_proj
      flag = 'proj_';
   else
       flag = '';
   end
   
   imwrite(im1Bary/255,[rep, 'Bary_', flag, filename{selec(1)}])
   imwrite(im2Bary/255,[rep, 'Bary_', flag, filename{selec(2)}])
   imwrite(im3Bary/255,[rep, 'Bary_', flag, filename{selec(3)}])
end


%% guided regularization

disp(' '), disp('Régularization ...'), disp(' ')
addpath('guided_filtering_code/')

N_gf = 10; % iterations
r = 4;
eps = 0.02^2;
q = zeros(N,M,3);

    
figure();
for i=1:K
    I = im{i};
    if     i==1, p = im1Bary;
    elseif i==2, p = im2Bary;
    elseif i==3, p = im3Bary;
    end
    I = I / 255; % guide = original image
    p = p / 255 - I; % image to be modified = OT color map

    for k=1:N_gf
        q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
        q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
        q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
        p = q;
    end
    subplot(1,K,i), imshow(q+I, [0, 1]);
    
    im_Bary_Reg{i} = q+I;
    
    if im_save && K==3
       rep = 'results/'
       if use_proj
          flag = 'proj_'; 
       else
           flag = '';
       end
   
   
       imwrite(im_Bary_Reg{i},[rep, 'Bary_Reg_', flag, filename{selec(i)}])
    end
end
clear I, clear P, clear q

%% Color transfer with multiple images (color grading)
if perform_grading
    
    disp(' '), disp('Color grading ...'), disp(' ')
    if K==3

        % parameters setting
        X1 = double(YCbCr{1});
        X2 = double(YCbCr{2});
        X3 = double(YCbCr{3});

        w = {                   [0,0,1]; ...
                            [1,0,3]; [0,1,3];  ...
                       [1,0,1]; [1,1,2]; [0,1,1]; ...
                   [3,0,1]; [2,1,1]; [1,2,1]; [0,3,1]; ...
                       [3,1,0]; [1,1,0]; [1,3,0]; [0,1,0]; ...
            }; % weights list
        options.method = 'grad'; % grad or stochastic
        Shift_Bary = zeros(2,numel(w));


        SW2_bary1 = [];
        tic, 
        figure
        for j=1:numel(w)
            % normalize weights
            [SW2_bary1{j},Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud(X1,{X1,X2,X3},w{j},options);

            im_b{j} = ycbcr2rgb( reshape(SW2_bary1{j}',[N M 3]) / 255 ) * 255;
            imshow(im_b{j}/255), drawnow
        end
        toc
    else
        return
    end


    % iterated guided regularization of the OT map
    N_gf = 5;
    r = 4;
    eps = 0.1^2;
    tic
    for i=1:numel(w)
        I = im{1};
        p = im_b{i}/255;

        I = I / 255; % guide = original image
        p = p - I; % image to be modified = OT color map

        for k=1:N_gf
            q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
            q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
            q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
            p = q;
        end

        im_b_r{i} = q + I;
    end
    toc
    clear I, clear P, clear q

    % creation d'une grande figure
    Z = zeros(N,M,3) + 255;
    Z2 = Z(:,1:round(M/2),:);
    Pic = cat(2,                   Z, Z, im_b_r{1}, Z, Z);
    Pic = cat(1, Pic, cat(2, Z2, Z, im_b_r{2}, im_b_r{3}, Z, Z2));
    Pic = cat(1, Pic, cat(2, Z, im_b_r{4}, im_b_r{5}, im_b_r{6}, Z));
    Pic = cat(1, Pic, cat(2, Z2, im_b_r{7}, im_b_r{8}, im_b_r{9}, im_b_r{10}, Z2));
    Pic = cat(1, Pic, cat(2, im{1}/255, im_b_r{11}, im_b_r{12}, im_b_r{13}, im_b_r{14}));
    figure,  imshow(Pic)

    if im_save
        imwrite(Pic,[rep, 'Proj_' filename{selec(1)}])
    end

end
    
return
