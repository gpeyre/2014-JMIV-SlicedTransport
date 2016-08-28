% Benchmark entre le code matlab et le code mex pour Sliced Wasserstein
% Projection

clear, close all,
% clc

Dim = 2; % dimension
n = 1e3; % number of points
d = 1e2; % number of directions

X0 = randn(Dim,n);
D = eye(Dim); D(1,1) = 2; D(2,2) = 1/2; 
Y = D*X0 + 3;

% test sur un cas critique
if 0 
    X0 = [1,0;0,1]
    Y = [1,0;1,0] + rand(2)*1e-2
end

w = [0 1]; w = w/sum(w); % weights for projection
t = (1:d)/d*pi; % uniform sampling or orientations
options.ndir = d;
options.base = [cos(t); sin(t)]; % directions set
options.step = 1; % descent step (maximum)
options.hessian = 0; % use hessian normalization
options.method = 'grad'; % 'grad' or 'stochastic'
options.nsubdir = d;
options.niter = 1e2;
options.display = 0;
options.eps = -1; % empeche l'ago de s'arreter à convergence pour une comparaison à nombre d'iteration fixé

tic, [X,E] = Sliced_Wasserstein_Projection_MeX(X0,Y, ...
                  options.step,options.base,options.niter, ...
                  strcmp(options.method, 'stochastic'),options.hessian,options.display); toc

figure, plot(X0(1,:),X0(2,:),'ok'), hold on, 
plot(Y(1,:),Y(2,:),'om'),
plot(X(1,:),X(2,:),'xb')
drawnow
%legend('Mex','Matlab')

tic, [X2,Energ2,E2] = Sliced_Wasserstein_Barycenter_PointCloud(X0,{X0,Y},w,options); toc

plot(X2(1,:),X2(2,:),'+r')

% energy monitoring
figure, semilogy(E,'b'), hold on, semilogy(Energ2,'r')
legend('Mex','Matlab')
