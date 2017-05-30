function [fout,gout,jac,r] = predictCSM_objective(kAlphaPhi, X, Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[fout,gout,jac,r]=f(x)
model = predictCSM(kAlphaPhi,X);
r = Y-model;
fout = .5*(r'*r);
if nargout>1
    x = X(:,1); y = X(:,2); w = X(:,3); 
    k = kAlphaPhi(1); alpha = kAlphaPhi(2); Phi = kAlphaPhi(3);
    tmp_in = -k*cos(alpha)*x-k*sin(alpha)*y+Phi;
    tmp_out = w.*[-cos(alpha)*x-sin(alpha)*y,k*sin(alpha)*x-k*cos(alpha)*y,ones(size(x,1),1)];
    jac_upper = tmp_out.*sin(tmp_in);
    jac_lower = -tmp_out.*cos(tmp_in);
    jac = [jac_upper;jac_lower];
    gout = jac'*r;
end
end

