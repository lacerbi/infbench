% FIGURE S? for VBMC 2020 paper. Plot g-and-k model results.

figname = {'vbmc2020_figSgandk'};
mypath = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
cd(mypath);

fontsize = 24;
axesfontsize = 20;
plotdim_flag = true;

options = [];
algos_list = [];

options.BestOutOf = 1;
options.NumZero = 1e-3;
options.Method = 'IR';
options.FontSize = fontsize;
options.AxesFontSize = axesfontsize;
options.ErrorBar = 1;
options.PlotType = 'line'; options.BootStrap = 0; options.Quantiles = [NaN,0.9];
options.PlotType = 'errorbar'; options.BootStrap = 1e5; options.Quantiles = [0.025,0.975];
algos_list{1} = {'wsabiplus@ldet','vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
algos_list{2} = {'parallelgp@v3','vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
algos_list{3} = {'parallelgp@v3','wsabiplus@ldet','vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
% algos = {'parallelgp@v3','vbmc@renewdefvarimiqrpluswup5fast','parallelgp'};

if 0
    algos_list = [];
    algos_list{1} = {'parallelgp@v3','vbmc@renewdefvarimiqrpluswup5fast','vbmc@viqrnogplate'};
    algos_list{1} = {'wsabiplus@ldet','parallelgp@v3'};
    options.BootStrap = 1e4;
end

probset = 'vbmc20';

nrows = 1;
ncols = 4;
panel_sub{1} = [1 2 3 4];
panel_legend = 4;
panel_empty = [];

labels = [];
labels{1} = 'A';
labels{2} = 'B';
labels{3} = 'C';

textpos = [0.75,0.9];

% Create figure, divide into panels
close all;

gutter = [.1, .12]; % H V
margins = [.1 .02 .2 .1]; % L R B T
hpanel = plotify(nrows,ncols,'gutter',gutter,'margins',margins,'labels',labels,'fontsize',fontsize);
for iPanel = panel_sub{1}(:)'; name{iPanel} = []; end

for iSub = 1:3    
    figure(1);
    
    algos = algos_list{min(iSub,numel(algos_list))};    
    options.DisplayLegend = false;
    
    switch iSub
        % case 1; options.Metric = 'costs'; options.NumZero = 0.1; options.YlimMax = 100;
        case 1; options.Metric = 'lnZ'; options.NumZero = 0.1; options.YlimMax = 1e3;
        case 2; options.Metric = 'mmtv'; options.NumZero = 0; options.YlimMax = 1;
        case 3; options.Metric = 'gsKL'; options.NumZero = 0.01; options.YlimMax = 1e4;
    end
    
    switch options.Metric
        case 'lnZ'; ystring = 'LML loss';
        case 'mmtv'; ystring = 'MMTV';
        case 'gsKL'; ystring = 'gsKL';
    end
    
    % Price (2018)
    dims = {'D1'};  noise = []; probs = {'price2018'};
    idx = panel_sub{1}(iSub);
    options.AxesHandles = hpanel(idx);
    options_temp = options;
    if strcmp(options.Metric,'lnZ'); options_temp.NumZero = 0.01; end
    if plotdim_flag && iSub == 1
        options_temp.PlotText = {textpos,'D = 4'};
    else
        options_temp.PlotText = [];
    end
    if iSub == 3
        options_temp.AxesHandles = [options_temp.AxesHandles,hpanel(panel_legend)];
        options_temp.DisplayLegend = true;
    end
    infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options_temp);    
    for iPanel = panel_sub{1}(iSub); axes(hpanel(iPanel))
        xlabel('Function evaluations','FontSize',fontsize'); 
        ylabel(ystring,'FontSize',fontsize); 
        set(gca,'FontSize',axesfontsize);
    end
end

%    for iPanel = panel_sub{iSub}(1:4); axes(hpanel(iPanel)); xlabel(''); end
name{1} = 'g-and-k';
 for iPanel = panel_sub{1}(:)'
     axes(hpanel(panel_sub{1}(iPanel)));
     if ~isempty(name{iPanel})
         title(name{iPanel},'FontSize',fontsize);
         set(gca,'FontSize',axesfontsize);
     end
 end

for iPanel = panel_empty
    axes(hpanel(iPanel));
    cla(gca,'reset');
    axis off;
end


% Save figure
pos = [1 41 1920 500];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;
print([mypath filesep() figname{1} '.pdf'],'-dpdf','-bestfit');

