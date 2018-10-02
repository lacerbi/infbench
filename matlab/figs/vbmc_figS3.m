% FIGURE S3 for VBMC paper. Plot example posterior of neuronal model.

names = {'x_1','x_2','x_3','x_4','x_5','x_6','x_7'};
fontsize = 18;

% Cornerplot of true posterior
id = 7;
infprob = infbench_goris2015([],id);
probstruct = infprob_init('vbmc18','goris2015',id,[],1,[]);

folder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data\goris2015mcmc';
cd(folder);
X = mergesamples(['goris2015_mcmc_n' num2str(id) '*'],[],0,0);
Xorig = warpvars(X,'inv',infprob.Data.trinfo);
bounds = [min(Xorig); max(Xorig)];
cornerplot(Xorig,names,[],bounds);
h = axes(gcf,'Position',[0 0 1 1]);
set(h,'Color','none','box','off','XTick',[],'YTick',[],'Units','normalized','Xcolor','none','Ycolor','none');
text(0.9,0.9,'True','FontSize',fontsize,'HorizontalAlignment','right');
title('Neuronal model (V2 neuron)','FontSize',fontsize);

pos = [100 200 980 780];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;

mypath = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
figname = 'vbmc_figS3a';
saveas(gcf,[mypath filesep() figname '.pdf']);

folder = ['C:\Users\Luigi\Dropbox\Postdoc\VBMC\data\vbmc18@goris2015\S' num2str(id)];
cd(folder);

close all;

% Cornerplot of variational posterior
n = 1;
data = load(['vbmc@base@' num2str(n) '.mat']);
stats = data.history{1}.Output.stats;
idx = find(stats.stable == 1,1);
vp = stats.vp(idx);

Xrnd = warpvars(vbmc_rnd(1e6,vp),'inv',probstruct.trinfo);
Xrnd = warpvars(Xrnd,'inv',infprob.Data.trinfo);
cornerplot(Xrnd,names,[],bounds);
h = axes(gcf,'Position',[0 0 1 1]);
set(h,'Color','none','box','off','XTick',[],'YTick',[],'Units','normalized','Xcolor','none','Ycolor','none');
text(0.9,0.9,['VBMC (iteration ' num2str(idx) ')'],'FontSize',fontsize,'HorizontalAlignment','right');

pos = [100 200 980 780];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;

figname = 'vbmc_figS3b';
saveas(gcf,[mypath filesep() figname '.pdf']);
