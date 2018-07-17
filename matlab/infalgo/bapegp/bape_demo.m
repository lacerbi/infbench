%BAPE_DEMO Run demo of BAPE on shallow Rosenbrock function.
% Reference:
%   1. Kandasamy, K., Schneider, J., & Póczos, B. (2017). Query efficient 
%      posterior estimation in scientific experiments via Bayesian active 
%      learning. Artificial Intelligence, 243, 45-56.

options = [];
options.Algorithm = 'bape'; % Select BAPE algorithm
options.Plot = 1;           % Plot posterior and search points each iteration
options.Meanfun = 'const';  % Constant GP mean function
% options.Meanfun = 'negquad';  % Negative quadratic GP mean function

fun = @rosenbrock_test;     % This is a shallow Rosenbrock function
x0 = [0 0];                 % Starting point
LB = [-5 -5];               % Hard bounds - this makes the problem easier
UB = [5 5];
PLB = [-1 -1];              % Plausible bounds identify initial region
PUB = [1 1];

[X,y,exitflag,output] = bapegp(fun,x0,LB,UB,PLB,PUB,options);

% Now try without the bounds...
[X,y,exitflag,output] = bapegp(fun,x0,[],[],PLB,PUB,options);