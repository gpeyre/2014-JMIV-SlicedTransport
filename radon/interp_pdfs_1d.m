%author: Nicolas Bonneel

function [ y ] = interp_pdfs_1d( varargin )
% the first N arguments are the pdfs, the last N arguments their weights.
% usage: interp_pdfs4(pdf1, pdf2, pdf3,..., weight1, weight2, weight3, ...)
% with sum(weights_i)=1, and pdf_i a vector

n_pdfs = length(varargin)/2;

interpinvCDF = zeros(size(varargin{1}));
for i=1:n_pdfs
 cdf_i = cumsum(varargin{i});    
 invcdf_i = inverseFunc(cdf_i, 1.0);    
 interpinvCDF = interpinvCDF+varargin{n_pdfs+i}*invcdf_i;
end

interpCDF = inverseFunc(interpinvCDF, length(interpinvCDF) );
y = [0; diff(interpCDF)];

end

