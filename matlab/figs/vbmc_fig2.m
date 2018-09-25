% FIGURE 2 for revised VBMC paper. Plot benchmark on synthetic likelihoods.

options.BestOutOf = 1;
options.NumZero = 1e-4;
options.Method = 'IR';
options.ErrorBar = 1;
options.BootStrap = 1e4;
options.SampleFrequency = NaN;
plots = {'lnZ','gsKL'};

algos = {'vbmc','vbmc@acqusreg','wsabi','wsabi@mm','bbq','bmc','agp','bape@nqreg','smc','ais'};

% algos = {'wsabi','wsabi@mm','bbq','bmc','agp','bape@negquad','smc','ais','vbmc','vbmc@acqusreg'};
dims = {'2D','6D','10D'};
% dims = {'2D','4D','6D','8D','10D'};
noise = [];

n = 1;
probset = 'vbmc18';
probs = {'lumpy','studentt','cigar'};

figname = {'vbmc_fig2a','vbmc_fig2b'};
mypath = fileparts(mfilename('fullpath'));
mypath = '.';

for iPlot = 1:numel(plots)
    options.PlotType = plots{iPlot};
    options.DisplayLegend = iPlot == numel(plots);
    figure(iPlot);
    infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);
    pos = [50,50,900,750];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    saveas(gcf,[mypath filesep() figname{iPlot} '.pdf']);
end