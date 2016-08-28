function [D, F, perm] = OT_mosek(C,options)

% OT_mosek - solve the Optimal Transport Problem for cost matrix C
%                       (X and Y, with n bins, may be multi-dimensionnal)
%
% Use : 
%       [Dist, F, perm] = OT_mosek(C)
%
% Solve :
%       min <C,F> where C_ij = cost(X_i - Y_j) from point clouds X_i and  Y_j
% subject to the following constraints on the optimal flow F:
%       F * ones(n,1) = ones(n,1)
%       F'* ones(n,1) = ones(n,1)
%       F >= 0
%
% Total cost is Dist = <C,F>  
% Flow is s.t Y = F*X 
% Perm is s.t Y = X(perm)
%
%   Copyright (c) 2013 G. Peyre, J. Rabin
%
% % Test :
% n = 1e2;
% X = rand(2,n); Y = rand(2,n)*3+1;
% X2 = repmat(sum(X.^2,1)', [1 n]);
% Y2 = repmat(sum(Y.^2,1) , [n 1]); 
% costMat = X2 + Y2 -2*X'*Y; % L2^2 distance : figure, imagesc(costMat)
%
% tic, [D,F,perm] = OT_mosek(costMat); tic
% D
% P = Y*F';
% figure, plot(X(1,:),X(2,:),'bo'), hold on,
% plot(Y(1,:),Y(2,:),'ro'), legend('X','Y')
% line([X(1,:);P(1,:)],[X(2,:);P(2,:)]), % assignments
% plot(Y(1,perm),Y(2,perm),'m+'), 

load_mosek();

options.null = 0;

N = size(C,1); % number of points
if (size(C,2) ~= N)
    disp('error : input cost matrix should be square')
    return 
end

%% Design constraint matrix

%  Amin <= A * F(:) <= Amax
A = [   ComputeR(ones(N,1)); ... % fait la somme des coefficients de chaque ligne, soit F*ones(n,1)
        ComputeL(ones(1,N)); ... % fait la somme des coefficients de chaque colonne, soit ones(1,n)*F
    ];
     
Amax = [   
           ones(N,1); ...
           ones(N,1); ...
       ]; 

Amin = [   
           ones(N,1); ...
           ones(N,1); ...
       ];
   
% Fmin <= F <= Fmax
Fmin = [ sparse(N*N,1) ]; % F >= 0
Fmax = []; % F <= 1 : inutile d'ajouter cette condition


%%
% Setup Mosek variables.

prob.c = C(:);
prob.a = A;
prob.blc = Amin;
prob.buc = Amax;
prob.blx = Fmin;
prob.bux = Fmax;

%%
% Set parameters.

param = [];
% max number of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 100);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', 1e-12);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', 1e-12);
% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);

% Perform the optimization.
tic, [r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);toc
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
err.niter = res.info.MSK_IINF_INTPNT_ITER;
sol = res.sol.itr;
w   = sol.xx;
F = reshape( w, [N N] ); % flow from X to Y

D = C.*F; D = sum(D(:));

[val,perm] = max(F,[],2);


