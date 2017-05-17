function [lmpass,lmgradnorm,lmtime,nlpass,nlgradnorm,nltime] = ...
    test_lmffd_onmodel(X,Y,modelfun,gradfun,beta0,options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Fit the parameters
tic;
[beta,resid,jacob,histout,costdata] = ...
    levmarfit_fd(X, Y, modelfun, beta0, options);
lmtime = toc;

try
    tic;
    [betanl,resid,jacob] = ...
        nlinfit(X,Y,modelfun,beta0,options);
    nltime = toc;
catch
    nltime = NaN;
end

% Ensure analytic gradient has been reduced below the tolerance
lmgradnorm = norm(gradfun(beta));
lmpass = norm(gradfun(beta))<options.TolX*10^3;

try
    nlgradnorm = norm(gradfun(betanl));
    nlpass = nlgradnorm<options.TolX*10^3;
catch
    nlgradnorm = NaN; nlpass = 0;
end
end

