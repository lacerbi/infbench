% FIGURE S6 for VBMC paper. Plot control on Goris et al. (2015).

options.BestOutOf = 1;
options.NumZero = 1e-2;
options.Method = 'IR';
options.ErrorBar = 1;
plots = {'lnZ','gsKL'};
options.BootStrap = 1e5;
%options.Quantiles = [0.25,0.75];
options.SampleFrequency = NaN;

algos = {'wsabi','vbmc','vbmc@control'};
dims = {'S8','S7'};
noise = [];

n = 1;
probset = 'vbmc18';
probs = {'goris2015'};

figname = {'vbmc_figS6a','vbmc_figS6b'};
mypath = fileparts(mfilename('fullpath'));
mypath = '.';

YlimMax = [1e4,1e6];

for iPlot = 1:numel(plots)
    options.PlotType = plots{iPlot};
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