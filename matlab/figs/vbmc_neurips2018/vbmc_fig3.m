% FIGURE 3 for revised VBMC paper. Plot benchmark on Goris et al. (2015).

options.BestOutOf = 1;
options.NumZero = 1e-2;
options.Method = 'IR';
options.ErrorBar = 1;
options.BootStrap = 1e5;
options.SampleFrequency = NaN;
plots = {'lnZ','gsKL'};

algos = {'smc','ais','bmc','wsabi@base2','wsabi@mm2','bbq','agp','bape@nqreg','vbmc@acqusreg','vbmc@acqproreg'};
% We also ran 'bbq@marginal' (BBQ*), with no much difference from BBQ. 
% For WSABI here we increased the duration of hyperparameter optimization,
% which ends up improving the exploration in active sampling.

dims = {'S8','S7'};
noise = [];

n = 1;
probset = 'vbmc18';
probs = {'goris2015'};

%algos = {'vbmc@acqproponly','bape'};

figname = {'vbmc_fig3a','vbmc_fig3b'};
mypath = fileparts(mfilename('fullpath'));
mypath = '.';

YlimMax = [1e4,1e6];

for iPlot = 1:numel(plots)
    options.Metric = plots{iPlot};
    options.YlimMax = YlimMax(iPlot);
    options.DisplayLegend = iPlot == numel(plots);
    figure(iPlot);
    infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);
    pos = [20,20,900,450];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    saveas(gcf,[mypath filesep() figname{iPlot} '.pdf']);
end