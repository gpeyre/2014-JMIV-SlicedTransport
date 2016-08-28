% Matlab Implementation of Error diffusion using Floyd's and Steinberg's...
% filter weights.
a=(imread('cameraman.tif')); % rgb2gray
s= size(a);
% error=ones(s(1),s(2));
for i = 1 : s(1)
       for j = 1 : s(2)
           if(a(i,j) < 127 )
               b(i,j)=0;
           else 
               b(i,j) = 255;
           end;
           qerror = a(i,j) - b(i,j);
       if(j < s(2))
           a(i,j+1) =  ((7/16 *qerror)+a(i,j+1)); 
       end;
       if(i<s(1) && j > 1)
           a(i+1,j-1) = a(i+1,j-1) + (3/16 *qerror);
           a(i+1,j) = a(i+1,j) + (5/16 *qerror);
       end;
       if(j<s(2) && i<s(1))
           a(i+1,j+1) = a(i+1,j+1) + (1/16 *qerror);
       end;
       end;
end;
  
% subplot(2,1,1);
% imshow(a);
% title('Original Gray Scale Image');
% subplot(2,1,2);
imshow(b);
title('Praveen s Half toning Creation');
