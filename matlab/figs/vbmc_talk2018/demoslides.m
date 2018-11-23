% FIGURE 1 for revised VBMC paper. Demo of VBMC on banana function.

if ~exist('stats','var'); stats = []; end
plotbnd = [-3,-2; 3,6];
close all;

rng(0);

% Define prior
sigma2 = 9;     % Prior variance
priorfun = @(x) - 0.5*sum(x.^2,2)/sigma2 - 0.5*size(x,2)*log(2*pi*sigma2);

stats = vbmc_film2d(@(x) rosenbrock_test(x) + priorfun(x),stats,plotbnd);
