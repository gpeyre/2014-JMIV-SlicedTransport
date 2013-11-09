% author: Nicolas Bonneel
% 2d Radon Wasserstein Barycenter
% usage: barycenter = interp_pdfs(1, pdf0, pdf1, pdf2, ...., weight0, weight1, weight2, ... )

function [ y ] = interp_pdfs(fastslant, varargin )

addpath('../fast-slant-stack/SlantStackP');

N = size(varargin{1},1)
N_imgs = length(varargin)/2;

radon_pdfs = cell(N_imgs,1);
mini = cell(N_imgs,1);
angles = linspace(0,179,180); % if fastslant==0 only.

% normalizes
for i=1:N_imgs
    radon_pdfs{i}=varargin{i}(:,:)/sum(sum(varargin{i}(:,:)));
    
    if (fastslant==1)
        radon_pdfs{i} = max(0,real(FastSlantStack(radon_pdfs{i}(:,:))))*N*2/sqrt(2);;
    else
        radon_pdfs{i} = max(0,real(radon(radon_pdfs{i}(:,:),angles)));
    end
        
end

radonInterp = interp_pdfs_Nd_radon(radon_pdfs{:}, varargin{N_imgs+1:end});
if (fastslant==1)
    y = real(Inv_FastSlantStack(radonInterp))*sqrt(2)/(N*2);
else
    y = real(iradon(radonInterp,angles));
end