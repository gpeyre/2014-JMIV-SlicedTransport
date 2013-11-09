%% implementation of the Radon barycenter
%% author: Nicolas Bonneel, 2013

N = 128;         % size of pdfs
N_imgs = 3;      % number of pdfs
N_interp = 201;  % number of time steps
use_fast_slant_stack = 1; % matlab's radon or not


% generates N_imgs  NxN 2D pdfs
Imgs = zeros(N,N,N_imgs);
for i=1:N
    for j=1:N
        Imgs(i,j,1) = exp(- ((i-2.3*N/3)^2 + (j-0.5*N/3)^2)/(2*(N/13)^2))/(N/13) + exp(- ((i-2.2*N/3)^2 + (j-1.2*N/3)^2 + (i-2.2*N/3)*(j-1.2*N/3)*0.5)/(2*(N/10)^2))/(N/10);
        Imgs(i,j,2) = exp(- ((i-0.7*N/3)^2*4 + (j-2*N/3)^2)/(2*(N/5)^2))/(N/5) + exp(- ((i-N/2)^2 + (j-3*N/4)^2)/(2*(N/9)^2))/(N/9) + exp(- ((i-0.8*N/3)^2*2 + (j-3*N/4)^2)/(2*(N/6)^2))/(N/6);
        Imgs(i,j,3) = exp(- ((i-1*N/3)^2 + (i-1*N/3)*(j-1.3*N/3)*4 + (j-1.3*N/3)^2*6)/(2*(N/7)^2))/(N/7) + exp(- ((i-1*N/3)^2 + 0.6*(i-1*N/3)*(j-1*N/4) + (j-1*N/4)^2*3)/(2*(N/8)^2))/(N/8)  + exp(- ((i-0.6*N/3)^2*3 + 1.3*(i-0.6*N/3)*(j-1*N/5) + (j-1*N/5)^2)/(2*(N/8)^2))/(N/8) + exp(- ((i-1*N/4)^2*2 + (j-1*N/5)^2)/(2*(N/6)^2))/(N/6);
        %Imgs(i,j,4) = exp(- ((i-2.3*N/3)^2 + (j-2.3*N/3)^2*6)/(2*(N/7)^2))/(N/7) + exp(- ((i-2*N/3)^2 + (j-3*N/4)^2*3)/(2*(N/8)^2))/(N/8)  + exp(- ((i-2.6*N/3)^2*3 + (j-3*N/5)^2)/(2*(N/8)^2))/(N/8) + exp(- ((i-2*N/4)^2*2 + (j-3*N/5)^2)/(2*(N/6)^2))/(N/6);
    end;
end;

% normalizes
for i=1:N_imgs
    Imgs(:,:,i)=Imgs(:,:,i)/sum(sum(Imgs(:,:,i)));
end

% display the first one
cte = 100./max(max(Imgs(:,:,1)));
figure(1);
image(Imgs(:,:,1)*cte);
colormap('gray')

% precompute Radon transforms
tic;
radon_pdfs = cell(N_imgs,1);
if (use_fast_slant_stack)
    addpath('../fast-slant-stack/SlantStackP');
    for i=1:N_imgs
        radon_pdfs{i} = max(0, real(FastSlantStack(Imgs(:,:,i))))*N*2/sqrt(2);
    end;
else
    for i=1:N_imgs
        radon_pdfs{i} = radon(Imgs(:,:,i), 0:179);
    end
end
toc;


% generates a path from pdf1 to pdf2, then to pdf3, then to their
% barycenter
tic;
result = zeros( N, N, N_interp);
% matlabpool(8);
parfor i=1:N_interp
    
    t = (i-1)/(N_interp-1);    
    if (t<0.25)
        a = t/0.25;
        b = 1-a;
    else
        if (t<0.5)
            a = 1-(t-0.25)/0.25;
            b=0;
        else
            if (t<0.75)
                a = 0;
                b = (t-0.5)/0.25;
            else
                a = (t-0.75)/0.25 * 1/3.;
                b = 1. - (t-0.75)/0.25 * 2/3.;
            end
        end;
    end;
    
    radonInterp = interp_pdfs_Nd_radon(radon_pdfs{:}, a, b, 1.-a-b);
    
    if (use_fast_slant_stack)
        result(:,:,i) = real(Inv_FastSlantStack(radonInterp))*sqrt(2)/(N*2);
    else
        result(:,:,i) = iradon(radonInterp, 0:179, 'linear', 'Hann', 1., N);
    end;
    
end;
matlabpool close;
toc

% display result
for i=1:N_interp
    figure(1);
    image(squeeze(result(:,:, i)*cte));
    colormap('gray')
    anim(i) = getframe;
end
movie(anim,1,10);
