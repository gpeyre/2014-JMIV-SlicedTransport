%author: Nicolas Bonneel

function  y = interp_pdfs_Nd_radon( varargin)  
% input: radon transforms of N 2d pdfs, and N weights
n_pdfs = length(varargin)/2;
size_pdf = size(varargin{1});
slices = cell(n_pdfs,1);

y = zeros(size_pdf);
for i=1:size_pdf(2)  % for each slice
    for j=1:n_pdfs
        slices{j}=varargin{j}(:,i);  
    end
    y(:,i) = interp_pdfs_1d( slices{:}, varargin{(n_pdfs+1):end}); % interpolates the slice
end

end

