function [x] = invcovmv(cenum, f, reg_enum, reg_cor, reg_sig)

% solves Cm x = f, by solving the conjugate graident normal equations
% (CGNE), Cm^*Cm x = Cm^*f.  Cm^* = Cm, acts as
% regularization/preconditioner.  Because it's symmetric we can use the
% same mat-vec multiplying routine when computed the transpose.

% Cm is the covariance matrix stored as an enumerated image linking
% homogenous regions

% f is the imput vector stored as an image
% x is the output vector stored as an image

% Regularization of Inverse Problems By Heinz W. Engl, Martin Hanke, Andreas Neubauer

itmax = 100;
eps = 0.001;

it = 1;
x = zeros(size(f));

vecAx = covmv(cenum, x, reg_enum, reg_cor, reg_sig);
r = f - vecAx;

vecAx = covmv(cenum, r, reg_enum, reg_cor, reg_sig);
s = vecAx;
p = s;

delt_new = real(sum(sum(conj(s).*s)));
delt_o = delt_new;

while it < itmax && delt_new > (eps^2)*delt_o
    
    q = covmv(cenum, p, reg_enum, reg_cor, reg_sig);
    alpha = real(delt_new)/real(sum(sum(conj(q).*q)));
    x = x + alpha*p;
    r = r - alpha*q;
    s = covmv(cenum, r, reg_enum, reg_cor, reg_sig);
    delt_old = delt_new;
    delt_new = real(sum(sum(conj(s).*s)));
    beta = real(delt_new)/real(delt_old);
    p = s + beta*p;
    
%      imagesc(x),colorbar,title(num2str(it))
%      disp([delt_new (eps^2)*delt_o])
%      pause
end


% d = r;
% delt_new = real(sum(sum(conj(r).*r)));
% delt_o = delt_new;
% 
% while it < itmax && delt_new > (eps^2)*delt_o
%     
%     q = covmv(cenum, d, reg_enum, reg_cor, reg_sig);
%     alpha = real(delt_new)/real(sum(sum(conj(d).*q)));
%     x = x + alpha*d;
% %     if mod(it,50) == 0
% %         vecAx = covmv(cenum, x, reg_enum, reg_cor, reg_sig);
% %         r = f - vecAx;
% %     else
%         r = r-alpha*q;
% %    end
%     delt_old = delt_new;
%     delt_new = real(sum(sum(conj(r).*r)));
%     beta = real(delt_new)/real(delt_old);
%     d = r + beta*d;
%     it = it + 1;
%     imagesc(r),colorbar,title(num2str(it))
%     disp([delt_new (eps^2)*delt_o])
%     pause
% end




