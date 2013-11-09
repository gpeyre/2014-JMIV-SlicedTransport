%% Script to see the Hessian matrix structure
clear all, close all, clc

use_half_sphere = 1 % 
use_random_dir = 1

%% 2D case
d = 2;
n = 1e2;
if ~use_random_dir % regular directions
    t = (0:n-1)/n;
else % random directions uniformly sampled
    t = rand(1,n);
end
    
if use_half_sphere
    t = t*pi;
else
    t = t*2*pi;
end

o = [cos(t); sin(t)];
H = o*o'

% [P,D] = eig(H)


figure, 
if n<=1e2
    line([zeros(1,n);o(1,:)],[zeros(1,n);o(2,:)])
else
    plot(o(1,:),o(2,:),'o')
end
axis([-1 1 -1 1]), axis equal

%% 3D case
d = 3;
n = 1e2; % n^2 directions
if ~use_random_dir % regular directions
    t = (0:n-1)/n;
    u = linspace(-1,1,n);
    % area : sin(v)dtdv = dtdu = mesure uniforme sur la sphere
    u = acos(u);
    n2 = n*n;

    if use_half_sphere
        t = t*pi;
    else
        t = t*2*pi;
    end
    
    [t,u] = meshgrid(t,u);
    t = t(:)'; u = u(:)';

    o = [cos(t).*sin(u); sin(t).*sin(u); cos(u)];
    
    
else % random directions uniformly sampled
    n2 = n*n;
    o = randn(d,n2);
    if use_half_sphere
       o(3,:) = abs(o(3,:)); 
    end
    o = o./repmat(sum(o.^2,1).^(1/2),[d 1]); % normalization
end

H = o*o'

% [P,D] = eig(H)


figure, 
if n2<=1e2
    line([zeros(1,n2);o(1,:)],[zeros(1,n2);o(2,:)],[zeros(1,n2);o(3,:)])
else
    plot3(o(1,:),o(2,:),o(3,:),'o')
end
axis([-1 1 -1 1 -1 1]), axis equal

%% d-D case : uniform distribution on the hyper sphere S^{d-1}

d = 1e2;
n = 1e3;

o = randn(d,n);
Hb = o*o' ;% normaly distributed vector
Ib = eye(d)*n; % asymtotic limit

o = o./repmat(sum(o.^2,1).^(1/2),[d 1]); % normalization
H = o*o'; % uniform on the hypersphere
I = eye(d)*n/d; % asymtotic limit

figure, 
subplot(2,3,1), imagesc(H), colorbar, title('Hessian matrix H'), axis equal, axis off
subplot(2,3,2), imagesc(I), colorbar, title('n/d*eye(d)'), axis equal, axis off
subplot(2,3,3), imagesc(abs(H-I)/n*d), colorbar, title('abs diff / n*d'), axis equal, axis off
subplot(2,3,4), imagesc(Hb), colorbar, title('Hessian matrix H with unormalized vector'), axis equal, axis off
subplot(2,3,5), imagesc(Ib), colorbar, title('n*eye(d)'), axis equal, axis off
subplot(2,3,6), imagesc(abs(Hb-Ib)/n), colorbar, title('abs diff / n'), axis equal, axis off


