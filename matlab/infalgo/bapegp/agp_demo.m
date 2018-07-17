%AGP_DEMO Run demo of AGP on shallow Rosenbrock function (as in [1]).
% Reference:
%   1. Wang, H., & Li, J. (2017). Adaptive Gaussian process approximation 
%      for Bayesian inference with expensive likelihood functions. 
%      arXiv preprint arXiv:1703.09930. 

options = [];
options.Algorithm = 'agp';  % Select AGP algorithm
options.Plot = 1;           % Plot posterior and search points each iteration

fun = @rosenbrock_test;     % This is a shallow Rosenbrock function
x0 = [0 0];                 % Starting point
LB = [-5 -5];               % Hard bounds - this makes the problem easier
UB = [5 5];
PLB = [-1 -1];              % Plausible bounds identify initial region
PUB = [1 1];

[X,y,exitflag,output,vbmodel] = bapegp(fun,x0,LB,UB,PLB,PUB,options);

% Now try without the bounds...
[X,y,exitflag,output,vbmodel] = bapegp(fun,x0,[],[],PLB,PUB,options);