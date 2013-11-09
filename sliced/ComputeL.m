function L = ComputeL(U, q, use_slow)

% ComputeL - left matrix multiplication as a linear operator
%
%   L = ComputeL(U, q, use_slow);
%
% Implement 
%       reshape(L(U)*sigma, [p,q]) = U*reshape(sigma,[n,q])
% where
% 	[p,n] = size(U);
%   [n,q] = size(Sigma);
% 
%   L has size (p*q,n*q)
%
%   Copyright (c) 2012 Gabriel Peyre


[p,n] = size(U);

if nargin<2
    q = n;
end
if nargin<3
    use_slow = 1;
end

% it is important to pre-alloc the size

if use_slow
    % slow code
    Q = length(find(U));
    L = spalloc(p*q,n*q,Q*q);
    for i=1:q
        L( 1+(i-1)*p:i*p, 1+(i-1)*n:i*n ) = U;
    end
else
    % fast code
    [i0,j0,v0] = find(U);
    [k,i] = meshgrid(0:q-1,i0(:)); i = i(:) + k(:)*p;
    [k,j] = meshgrid(0:q-1,j0(:)); j = j(:) + k(:)*n;
    v = repmat(v0(:), [q 1]);
    L = sparse(i,j,v(:), p*q, n*q);
end
