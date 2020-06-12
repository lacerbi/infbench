% FIGURE 2, 3 and 4 for VBMC 2020 paper. Plot benchmark on all test functions.

figname = {'vbmc2020_Sq90lml','vbmc2020_Sq90mmtv','vbmc2020_figSq90gskl'};
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
% plots = {'lnZ','gsKL','mmtv'};
options.FontSize = fontsize;
options.AxesFontSize = axesfontsize;
options.Metric = 'mmtv';
options.ErrorBar = 1;
options.PlotType = 'line'; options.BootStrap = 0; options.Quantiles = [NaN,0.9];
options.PlotType = 'errorbar'; options.BootStrap = 1e5; options.Quantiles = [0.025,0.975];
options.BaseQuantile = 0.9;
algos_list{1} = {'wsabiplus@ldet','vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
algos_list{2} = {'parallelgp@v3','vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
% algos = {'parallelgp@v3','vbmc@renewdefvarimiqrpluswup5fast','parallelgp'};

if 0
    algos_list = [];
    algos_list{1} = {'parallelgp@v3','vbmc@renewdefvarimiqrpluswup5fast','vbmc@viqrnogplate'};
    algos_list{1} = {'parallelgp@v3'};
    options.BootStrap = 0; options.Quantiles = [0.5,0.9];
end

probset = 'vbmc20';

nrows = 2;
ncols = 5;
panel_sub{1} = [1 2 3 4 5, 6 7 8 9 10];
panel_sub{2} = [1 2 3 4 5, 6 7 8 9 10];
panel_sub{3} = [1 2 3 4 5, 6 7 8 9 10];
panel_legend = 5;
panel_empty = [];

textpos = [0.75,0.9];

close all;
for iSub = 2
    figure(iSub);
    
    for iPanel = panel_sub{iSub}(:)'; name{iPanel} = []; end
    
    % Create figure, divide into panels
    gutter = [.06, .12]; % H V
    margins = [.07 .02 .1 .05]; % L R B T
    hpanel = plotify(nrows,ncols,struct('gutter',gutter,'margins',margins));
    algos = algos_list{min(iSub,numel(algos_list))};    
    
    options.DisplayLegend = false;
    
    switch iSub
        % case 1; options.Metric = 'costs'; options.NumZero = 0.1; options.YlimMax = 100;
        case 1; options.Metric = 'lnZ'; options.NumZero = 0.1; options.YlimMax = 1e3;
        case 2; options.Metric = 'mmtv'; options.NumZero = 0; options.YlimMax = 1;
        case 3; options.Metric = 'gsKL'; options.NumZero = 0.01; options.YlimMax = 1e4;
    end
    
    switch options.Metric
        case 'lnZ'; ystring = 'LML loss (90% quantile)';
        case 'mmtv'; ystring = 'MMTV (90% quantile)';
        case 'gsKL'; ystring = 'gsKL (90% quantile)';
    end
    
    % Wood (2010)
    dims = {'D1'};  noise = []; probs = {'wood2010'};
    idx = panel_sub{iSub}(1);
    options.AxesHandles = hpanel(idx);
    name{1} = 'Ricker';
    if plotdim_flag; options.PlotText = {textpos,'D = 3'}; else; options.PlotText = []; end
    data = infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);

    if 1
        if 0
            % Price (2018)
            dims = {'D1'};  noise = []; probs = {'price2018'};
            idx = panel_sub{iSub}(2);
            options.AxesHandles = hpanel(idx);
            options_temp = options;
            if strcmp(options.Metric,'lnZ'); options_temp.NumZero = 0.01; end
            name{2} = 'g-and-k';
            if plotdim_flag; options_temp.PlotText = {textpos,'D = 4'}; else; options_temp.PlotText = []; end
            infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options_temp);
        end
        
        % Krajbich et al. (2010)
        dims = {'S1','S2'};  noise = []; probs = {'krajbich2010'};
        idx = panel_sub{iSub}(2:3);        
        options.AxesHandles = hpanel(idx);
        name{2} = 'aDDM (S1)';
        name{3} = 'aDDM (S2)';
        if plotdim_flag; options.PlotText = {textpos,'D = 4'}; else; options.PlotText = []; end
        infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);
                
        % Acerbi (2012)
        dims = {'S1'};  noise = []; probs = {'acerbi2012'};
        idx = panel_sub{iSub}(4);
        options.AxesHandles = hpanel(idx);
        options_temp = options;
        name{4} = 'Timing';
        if plotdim_flag; options_temp.PlotText = {textpos,'D = 5'}; else; options_temp.PlotText = []; end
        infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options_temp);
                
        % Acerbi, Dokka, et al. (2018)
        dims = {'S1','S2'};  noise = []; probs = {'acerbidokka2018'};
        idx = panel_sub{iSub}(6:7);
        options.AxesHandles = hpanel(idx);
        name{6} = 'Multisensory (S1)';
        name{7} = 'Multisensory (S2)';
        if plotdim_flag; options.PlotText = {textpos,'D = 6'}; else; options.PlotText = []; end
        infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);

        % Akrami et al. (2018)
        dims = {'S1'};  noise = []; probs = {'akrami2018b'};
        idx = panel_sub{iSub}(10);
        options.AxesHandles = hpanel(idx);
        name{10} = 'Rodent';
        if plotdim_flag; options.PlotText = {textpos,'D = 9'}; else; options.PlotText = []; end
        infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options);

        % Goris et al. (2015)
        dims = {'S108','S107'};  noise = 'me'; probs = {'goris2015b'};
        idx = panel_sub{iSub}(8:9);
        options.AxesHandles = [hpanel(idx),hpanel(panel_legend)];
        options.DisplayLegend = true;
        name{8} = 'Neuronal (V1)';
        name{9} = 'Neuronal (V2)';
        if plotdim_flag; options.PlotText = {textpos,'D = 7'}; else; options.PlotText = []; end
        options_temp = options;
        if strcmp(options.Metric,'gsKL'); options_temp.YlimMax = 1e6; end        
        infbench_plot(probset,probs,dims,noise,algos,[],{'prob','subprob'},options_temp);
    end
    
    for iPanel = panel_sub{iSub}(1:4); axes(hpanel(iPanel)); xlabel(''); end
    for iPanel = panel_sub{iSub}(6:10); axes(hpanel(iPanel)); xlabel('Function evaluations','FontSize',fontsize'); end
    for iPanel = panel_sub{iSub}([1,6]); axes(hpanel(iPanel)); ylabel(ystring,'FontSize',fontsize); end
    for iPanel = panel_sub{iSub}(:)'
        axes(hpanel(panel_sub{iSub}(iPanel)));
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
    pos = [1 41 1920 963];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    print([mypath filesep() figname{iSub} '.pdf'],'-dpdf','-bestfit');
    
end


