function [vec] = covmv(cenum, im, reg_enum, reg_cor, reg_sig)

% computes the matrix vector product of the covariance matrix with a vector
% the covariance matrix is stored as an enumerated image relating homogeneous regions
% the input and output vectors (vec and im) are images stored as 2D or 3D


% cenum:    enumerated image of homogenous regions, e.g.
%             cenum = [1 1 1 1 1 4;
%                      1 1 2 2 1 4;
%                      1 1 2 2 1 1;
%                      1 1 1 1 3 3];
% im:       image 
% reg_enum: list of enumerations, e.g. reg_enum = [1 2 3 4];, with
%           corresponding correlations and sigmas below...
% reg_cor:  correlations of each region, e.g. cor_enum = [0.1 0.1 0.2 0.2];
% reg_sig:  standard deviation of pixels in each region, e.g. reg_sig = [10 10 10 2]


nreg = length(reg_enum);
vec = zeros(size(im));

for reg = 1:nreg,
    ind = find(cenum == reg_enum(reg));
    
    sig = reg_sig(reg);
    cor = reg_cor(reg);
    
    vec(ind) = sig^2*cor*sum(im(ind)) + sig^2*(1-cor)*im(ind);
end
