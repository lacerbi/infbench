% FIGURE S2 for revised VBMC paper. Plot full benchmark on synthetic likelihoods.

options.BestOutOf = 1;
options.NumZero = 1e-4;
YlimMax = [1e4,1e6];
options.Method = 'IR';
options.ErrorBar = 1;
options.BootStrap = 1e5;
options.SampleFrequency = NaN;
options.AdaptiveYlim = 1;
plots = {'lnZ','gsKL'};

%algos = {'smc','ais','bmc','wsabi','wsabi@mm','bbq','agp','bape@nqreg','vbmc@acqusreg','vbmc'};
algos = {'vbmc@acqproreg','smc','bmc','wsabi@base2','bbq','bape@nqreg'};

dims = {'2D','4D','6D','8D','10D'};
noise = [];

n = 1;
probset = 'vbmc18';
probs = {'lumpy','studentt','cigar'};

figname = {'vbmc_figS2a','vbmc_figS2b'};
mypath = fileparts(mfilename('fullpath'));
mypath = '.';

for iPlot = 1:numel(plots)
    options.Metric = plots{iPlot};
    options.YlimMax = YlimMax(iPlot);
    options.DisplayLegend = iPlot == numel(plots);
    figure(iPlot);
    infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);
    pos = [20,20,1400,1000];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    saveas(gcf,[mypath filesep() figname{iPlot} '.pdf']);
end