% FIGURE 1 for revised VBMC paper. Demo of VBMC on banana function.

if ~exist('stats','var'); stats = []; end
close all;

rng(0);

% Define prior
if 0
    sigma2 = 9;     % Prior variance
    priorfun = @(x) - 0.5*sum(x.^2,2)/sigma2 - 0.5*size(x,2)*log(2*pi*sigma2);
    plotbnd = [-3,-2; 3,6];
    jointfun = @(x) rosenbrock_test(x) + priorfun(x);
else
    M = 12;
    w = rand(M,1); w = w/sum(w);
    mu = rand(M,2)*5-2.5;
    jointfun = @(x) log(sum(w(:).*mvnpdf(x,mu,[1 0;0 1]),1));
    plotbnd = [-5,-5; 5,5];
end

stats = vbmc_film2d(jointfun,stats,plotbnd);
