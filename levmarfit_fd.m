function [beta,resid,jac,histout,costdata] = levmarfit_fd(X,Y,modelfun,beta0,options)
%
%LEVMARFIT_FD  nonlinear least-squares regression with Levenberg-Marquardt.
% [beta,resid,jac,histout,costdata] = LEVMARFIT_FD(X,Y,modelfun,beta0,options)
%
% Calling sequence based on the NLINFIT function for nonlinear regression.
% Repackages as a minimization problem for LEVMAR (C.T. Kelley)
%   with finite difference Jacobian.
%
% Input: X = predictor variables
%        Y = observations (column vector) with one row for each row of X
%        modelfun = modelfun(beta,X) is the model to be fit, where
%                     beta is a vector of parameters
%        beta0 = initial guess for optimum parameters
%        options = struct of options - must have TolX and MaxIter set
%
% Output: beta = solution
%         resid = final residual
%         jac = final Jacobian
%         histout = iteration history - each row of histout is      
%                [norm(grad), f, number of stepsize cuts, iteration count] 
%         costdata = [num f, num grad, num hess] (for levmar, num hess=0)
%

% Create objective function for LEVMAR
    f = @(beta) objective_lmffd(beta,X,Y,modelfun);
% Set tolerance and maximum iterations
    tol = options.TolX;
    maxit = options.MaxIter;
% Minimize with LEVMAR
    [beta,histout,costdata] = levmar(beta0,f,tol,maxit);
% Get Jacobian and residual with one more call to objective function
    [fout,gout,jac,resid] = f(beta);
end

function [fout,gout,jac,r] = objective_lmffd(beta,X,Y,modelfun)
    % calculate r(beta) and fout
    m = modelfun(beta,X);
    r = Y-m;
    fout = r'*r/2;
    if nargout>1
        % calculate r'(beta) (jac) and gout
        dbeta = eps^(1/3);
        jac = zeros(size(r,1),size(beta,1));
        h = dbeta*eye(size(beta,1));
        for i = 1:size(beta,1)
            mh = modelfun(beta+h(:,i),X);
            % Forward difference (r' = -m')
            jac(:,i) = (m-mh)/(dbeta);
        end
        gout = jac'*r;
    end
end