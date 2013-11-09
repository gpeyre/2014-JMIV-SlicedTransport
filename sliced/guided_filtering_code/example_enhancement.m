% example: detail enhancement
% figure 6 in our paper

close all;
tic

%I = double(imread('img_enhancement/tulips.bmp')) / 255;


rep = '/users/rabin/doc/matlab/IMAGE/'; % ls(rep)
name = [rep 'barbara_small.jpg' ] % 'kodak/kodim03.png' budgerigar_small.jpg lighthouse.png barbara_small.png cezanne Lena barbara_gray barbara Cameraman fingerprint house man
I = double(imread(name)) / 255;

%u = u(1:512, 1:512, :);

p = I;

r = 16;
eps = 0.1^2;

q = zeros(size(I));

q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);

I_enhanced = (I - q) * 5 + q;

toc

figure();
imshow([I, q, I_enhanced], [0, 1]);
