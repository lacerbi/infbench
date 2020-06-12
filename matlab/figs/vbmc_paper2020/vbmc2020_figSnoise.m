% FIGURE S? for VBMC 2020 paper. Plot noise-sensitivity results.

figname = {'vbmc2020_figSnoise'};
mypath = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
cd(mypath);

maxfiles = 50;

bigfontsize = 24;
fontsize = 18;
axesfontsize = 14;

nrows = 1;
ncols = 7;
panel_sub = [1 2 3 4 5 6 7];
panel_legend = 7;
panel_empty = [];

labels = [];
labels{1} = 'A';
labels{5} = 'B';

titlestring = [];
titlestring{1} = 'aDDM (S1)';
titlestring{5} = 'aDDM (S1)';
titlestring{2} = 'Timing';
titlestring{6} = 'Timing';
titlestring{3} = 'Neuronal (V1)';
titlestring{7} = 'Neuronal (V1)';

 
textpos = [0.75,0.9];

% Create figure, divide into panels
close all;

gutter = [.05, .12]; % H V
margins = [.1 .02 .2 .1]; % L R B T
hpanel = plotify(nrows,ncols,'gutter',gutter,'margins',margins,'labels',labels,'fontsize',bigfontsize);
for iPanel = panel_sub(:)'; name{iPanel} = []; end

noise_plot{1} = {'krajbich2010','S101',1};
noise_plot{5} = {'krajbich2010','S101',2};
noise_plot{2} = {'acerbi2012','S101',1};
noise_plot{6} = {'acerbi2012','S101',2};
noise_plot{3} = {'goris2015b','S108',1};
noise_plot{7} = {'goris2015b','S108',2};


% noise_plot{1} = {'krajbich2010','S102',1};
% noise_plot{5} = {'krajbich2010','S102',2};
% noise_plot{2} = {'akrami2018b','S101',1};
% noise_plot{6} = {'akrami2018b','S101',2};
% noise_plot{3} = {'goris2015b','S107',1};
% noise_plot{7} = {'goris2015b','S107',2};
% noise_plot{4} = [];
ylabel_idx = [1,5];
legend_idx = 7;

for iPanel = 1:numel(noise_plot)
    figure(1);
    axes(hpanel(iPanel));
    if isempty(noise_plot{iPanel})
        cla(gca,'reset');
        axis off;
        continue; 
    end
    
    prob = noise_plot{iPanel}{1};
    subprob = noise_plot{iPanel}{2};
    metric_idx = noise_plot{iPanel}{3};
    vbmc20_performance_noise(prob,subprob,[],metric_idx,maxfiles);
    if ~any(iPanel == ylabel_idx)
        xlabel(''); 
        ylabel(''); 
        % set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
    % set(gca,'XTick',0:1:5);
    if ~any(iPanel == legend_idx); legend off; end
    if ~isempty(titlestring{iPanel})
        title(titlestring{iPanel},'FontSize',fontsize);
    end
end
cd(mypath);

% Save figure
pos = [1 41 1920 500];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;
print([mypath filesep() figname{1} '.pdf'],'-dpdf','-bestfit');
