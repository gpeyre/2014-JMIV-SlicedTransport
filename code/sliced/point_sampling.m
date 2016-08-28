function X = point_sampling(im,Nb)
%
% Usage : 
%           X = point_sampling(im,N);
% Inputs :
%           im : input image used as a density
%           Nb : number of points sampled according to empirical probability distribution
% Outputs      
%           X  : point cloud coordinates [2 x Nb]
%
% Julien.rabin@unicaen.fr - 2013 (c)
%
% % For testing :
% clear, clc, close all
% Nb = 1e5;
% if 1
%     im = imread('cameraman.tif');
% elseif 1
%     n = 1e2;
%     im = zeros(2*n);
%     im(1:n,1:n) = 255;
%     im(n+1:2*n,n+1:2*n) = 255;
% else
%     im = zeros(1e2);
%     im(:,1:50) = 255;
% end
% figure(), imshow(im)
%
% % upsampling
% k = 2;
% I = imresize(im, k);
% tic, X = point_sampling(I,Nb); toc
% 
% figure, imshow(I)
% hold on, plot(X(1,:),X(2,:),'.'); axis off
%
% d = density_estimate(X,size(im,1),3);
% figure, imagesc(d), axis equal, axis off, colormap gray

%% 1) compute empirical distribution and repartition function
im = double(im)/255;

[n,m,d] = size(im);
if d>1
    disp('error : im should be a 2D array.')
    return 
end
N = n*m;

% disp('Normaliser image entre [0,1] ou non ?')
v = im(:); v = v / sum(v(:));

% sort values 
b = randn(N,1)*0e-9; % additive noise to have unique pixel values
[v,R] = sort(v+b,'ascend');
if 0 % Compute inverse ranking function (useless here)
    Id = (1:N); % identity
    Rinv = zeros(N,1);
    Rinv(R) = Id; % reverse rank index
    Rinv = Rinv(:);
end

% empirical repartition function
F = cumsum(v);
% F = (F-F(1))/(F(end)-F(1)); 
% figure, plot(F,'x-')



%% 3) Choose random pixel accordingly to the input PdF
r = rand(Nb,1); % uniform random values
r = sort(r,'ascend');

% pseudo inverse de F
Idx = zeros(Nb,1);
j=1;
for i = 1 : Nb
   while r(i) > F(j)
      j = j+1; 
   end
   Idx(i) = j;
end
% figure, plot(Idx,'x-')
% figure, plot(R(Idx),'x-')

[I,J] = ind2sub([n m],R(Idx));

% figure(Nfig), hold on, plot(J(:),I(:),'.'); axis off

% disp('ToDo : error diffusion scheme on repartition function')

X = [J(:)';I(:)'];

return
