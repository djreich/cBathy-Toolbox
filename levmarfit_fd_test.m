try
    options = statset('nlinfit');
catch
end
options.MaxIter = 100;
options.TolX = 1e-10;
options.Display = 'off';
% Produce the data
XX = cell(1,4);
XX{1} = (0:0.01:1)'; XX{2} = (0:1:100)'; 
XX{3} = (0:0.0001:1)'; XX{4} = (0:0.01:100)';
YY = cell(1,4);

%% First test case: 
% Fit Y= a sin( bX ) + b sin( aX )
% Exact solution is (a = 2 , b = 0.5)

% Function to fit
modelfun=@(beta,X)beta(1)*sin(beta(2)*X)+beta(2)*cos(beta(1)*X);
% Produce observations
for i=1:4; YY{i}=modelfun([2;0.5],XX{i}); end
% Track timings
timing = zeros(4,2);
% Test 4 cases
for i=1:4
    X = XX{i}; Y = YY{i};
    % Analytic gradient (of 1/2 r'r) for verification
    gradfun=@(beta)...
        [(modelfun(beta,X)-Y)'*(sin(beta(2)*X)-beta(2)*X.*sin(beta(1)*X));
        (modelfun(beta,X)-Y)'*(beta(1)*X.*cos(beta(2)*X)+cos(beta(1)*X))];
    beta0 = [2-1/X(end);0.5+1/X(end)];
    [lmpass,lmgradnorm,lmtime,nlpass,nlgradnorm,nltime] = ...
        test_lmffd_onmodel(X,Y,modelfun,gradfun,beta0,options);
    assert(lmpass,...
        sprintf('Levenberg-Marquardt failed to reduce gradient norm below tolerance in data case %d',i));
    timing(i,1)=lmtime; timing(i,2)=nltime;
end
disp('timings [levmarfit_fd,nlinfit]')
disp(timing);

%% Second test case: Rosenbrock's function
% Minimize f(x,y) = (1-x)^2 + 100(y-x^2)^2
% Exact solution is (x = 1, y = 1)

% Function to fit
modelfun=@(beta,X)(1-beta(1))^2+100*(beta(2)-beta(1)^2)^2;
X = 0; Y = 0;
gradfun=@(beta)...
    modelfun(beta,X)*[-2*(1-beta(1))-400*(beta(2)-beta(1)^2)*beta(1);
                       200*(beta(2)-beta(1)^2)];
beta0 = [1.01;0.99];
[lmpass,lmgradnorm,lmtime,nlpass,nlgradnorm,nltime] = ...
    test_lmffd_onmodel(X,Y,modelfun,gradfun,beta0,options);
assert(lmpass,...
    sprintf('Levenberg-Marquardt failed to reduce gradient norm below tolerance in data case %d',i));

%% Third test case: Regression Test on Single Bathy Run (fCombined)
%disable Statistics Toolbox
rmpath(strcat(matlabroot,'/toolbox/stats/stats/'));
stationStr = 'argus02a';
stackName = 'testStack102210Duck';
%get bathy struct without using Stat toolbox
bathy_ns = analyzeSingleBathyRunNotCIL(stackName, stationStr);
%re-enable Statistics Toolbox
addpath(strcat(matlabroot,'/toolbox/stats/stats/'));
%get bathy struct using Stat toolbox
try
    load('bathystruct-102210Duck.mat');
catch
    bathy = analyzeSingleBathyRunNotCIL(stackName, stationStr);
end
h = bathy.fCombined.h; hns = bathy_ns.fCombined.h;
% fCombined depth estimates have same NaN structure
assert(isequal(isnan(h),isnan(hns)),'NaN values appear in different locations in fCombined.h')
diff = abs(h(~isnan(h))-hns(~isnan(hns)))./(eps+abs(h(~isnan(h))))>1e-3;
% fCombined depth estimates differ by no more than 0.1 percent
assert(sum(diff)==0,sprintf('fCombined.h does not match in %d locations',sum(diff)));
hErr = bathy.fCombined.hErr; hnsErr = bathy_ns.fCombined.hErr;
% hErr has same NaN structure
assert(isequal(isnan(hErr),isnan(hnsErr)),'NaN values appear in different locations in fCombined.hErr');
diffErr = abs(hErr(~isnan(hErr))-hnsErr(~isnan(hnsErr)))./(eps+abs(hErr(~isnan(hErr))))>1e-3;
% hErr differ by no more than 0.1 percent
assert(sum(diffErr)==0,sprintf('fCombined.hErr does not match in %d places',sum(diffErr)));
J = bathy.fCombined.J; Jns = bathy_ns.fCombined.J;
% J has same NaN structure
assert(isequal(isnan(J),isnan(Jns)),'NaN values appear in different locations in fCombined.J');
diffJ = abs(J(~isnan(J))-Jns(~isnan(Jns)))./(eps+abs(J(~isnan(J))))>1e-3;
% J differs by no more than 0.1 percent
assert(sum(diffJ)==0,sprintf('fCombined.J does not match in %d locations',sum(diffJ)));