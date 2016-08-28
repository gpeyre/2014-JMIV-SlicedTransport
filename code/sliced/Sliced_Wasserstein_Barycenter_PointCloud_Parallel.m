function [X,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud_Parallel(X0,Y,w,options)

% Sliced Wasserstein Barycenter - multi-dimensional barycenter using sliced wasserstein energy descent
%
%   [X,Energ,E] = Sliced_Wasserstein_Barycenter(X0,Y,w,options);
%
% OUTPUTs
%
%   X is the Sliced Wasserstein barycenter of Y, minimizing the following energy:
%       = min_{s_j} 1/2 * Sum_o Sum_j Sum_i { w(j) * < X(i)-Y{j}(s_j(i)) , o >^2 }
%
%   Energ(k) is the slice wasserstein distance at iteration k : SW(X{k+1}, X{k})
%
%   E(k) is the descent amount at iteration k : ||X{k+1} - X{k}||_2
%   
% INPUTs
%
%  X0 : initialization of barycenter X (should be one of Y distributions)
%
%  Y : list of N-point-wise distributions in R^D (cell array of D-by-N matrices)
%
%  options:
%   options.niter is the number of iterations (100 by default)
%   options.ndir is the number of directions (ndir=dimension by default)
%   options.nsubdir is the number of directions if stochastic descent (nsubdir=ndir by default)
%   options.step is the step descend (step=1 by default).
%   options.eps sets the threshold on normalized displacement/descent amount E, in order to stop prematurely the algorithm
%
%   options.ortho_base=1 uses Gram-Schmidt algorithm for the first 'd' vectors of the base (=dimension) 
%   options.random_base=1 uses random base for each step (mandatory for stochastic gradient descent)
%           else if options.base = [d x ndir] matrix, it is used as fixed Bases
%
%   options.method is the descent strategy : 
%     - 'grad' : gradient descent (classical averaging batch descent which can be combined with Hessian 'normalization')
%     - 'stochastic' : stochastic gradient descent on a subset of directions
%
%   options.hessian=1 uses hessian matrix to scale gradient step
%
%   options.display=1 displays intermediate states (2D/3D only)
%   options.message=1 allows printing
%   options.pause : = 1 adds a pause, = t adds extra time of t seconds at each iteration, 
%
%
% EXAMPLE
%
%     n = 1e3;
%     x = linspace(0,1,n)*2*pi;
%     Y = { [cos(x); sin(x); ]; [cos(x)+4; sin(2*x); ]; [cos(2*x)+2; sin(3*x)+2; ]; };
%     nb = figure; plot(Y{1}(1,:), Y{1}(2,:), 'r'), hold on
%     plot(Y{2}(1,:), Y{2}(2,:), 'g'), plot(Y{3}(1,:), Y{3}(2,:), 'b')
%
%     X0 = Y{1};
%     w = [1 1 1];
%     d = 1e2;
%     t = (1:d)/d*pi;
%     options.ndir = d;
%     options.method = 'stochastic';
%     options.nsubdir = round(d/10);
%     options.step = 1;
%     options.display = 1;
%     options.base = [cos(t); sin(t)];
%     tic, [X,Energ,E] = Sliced_Wasserstein_Barycenter_PointCloud_Parallel(X0,Y,w,options); toc 
%     figure(nb), plot(X(1,:), X(2,:), '.k')
%
%
% Copyright, J. Rabin & G. Peyré , 2013.
% CNRS - Cérémade, Dauphine University, Paris, France
%
% See Also SLICED_WASSERSTEIN_PROJECTION, SLICED_WASSERSTEIN_KMEANS
%

disp('ToDo')
