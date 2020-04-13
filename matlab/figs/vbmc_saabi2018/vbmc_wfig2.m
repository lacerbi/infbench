% FIGURE 2 for VBMC workshop paper. Plot benchmark on Goris et al. (2015).

options.BestOutOf = 1;
options.NumZero = 1e-1;
options.Method = 'IR';
options.ErrorBar = 1;
options.BootStrap = 1e5;
options.SampleFrequency = NaN;
options.AdaptiveYlim = 1;
plots = {'lnZ','gsKL'};

algos = {'vbmc@acqusreg','vbmc@acqf2reg','vbmc@acqfregvlnn','vbmc@acqfregvsqrtn','vbmc@const','vbmc@se','vbmc@acqproreg'};
dims = {'S8','S7'};
noise = [];

n = 1;
probset = 'vbmc18';
probs = {'goris2015'};

%algos = {'vbmc@acqproponly','bape'};

figname = {'vbmc_wfig2a','vbmc_wfig2b'};
mypath = fileparts(mfilename('fullpath'));
mypath = '.';

YlimMax = [1e3,1e4];

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