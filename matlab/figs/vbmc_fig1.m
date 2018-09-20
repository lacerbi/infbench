% FIGURE 1 for revised VBMC paper. Demo of VBMC on banana function.

if ~exist('stats','var'); stats = []; end
plotbnd = [-3,-2; 3,6];
close all;
stats = vbmc_demo2d(@rosenbrock_test,stats,plotbnd);

pause

mypath = '.';
figname = 'vbmc_fig1';
saveas(gcf,[mypath filesep() figname '.pdf']);
