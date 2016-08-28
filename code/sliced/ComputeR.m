function R = ComputeR(V, q, use_slow)

% ComputeR - right matrix multiplication as a linear operator
%
%   R = ComputeR(V, q, use_slow);
%
% Implement 
%       reshape(R(V)*sigma, [q,p]) = reshape(sigma,[q,n])*V
% where
% 	[n,p] = size(V);
%   [q,n] = size(Sigma);
%   [q,p] = size(Sigma*V)
%
%   Copyright (c) 2012 Gabriel Peyre

[n,p] = size(V);

if nargin<2
    q = n;
end
if nargin<3
    use_slow = 1;
end

% Use 
%  Sigma*V = (V' * Sigma')'
%  L(V) = perm(q,p) o L(V') o perm(q,n)

L = ComputeL(V', q, use_slow);
R = perm_mat(p,q) * L * perm_mat(q,n);


% I = reshape(1:n*p, [n p])'; I = I(:);
% J = reshape(1:n*n, [n n])'; J = J(:);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = perm_mat(p,n)

% permute applied to (p,n) matrix
[b,a] = meshgrid(1:n,1:p); a = a(:); b = b(:);
P = sparse(b + (a-1)*n, a + (b-1)*p, ones(n*p,1));

end
