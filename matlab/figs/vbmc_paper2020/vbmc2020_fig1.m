% FIGURE 1 for VBMC 2020 paper. 
% Demo of different acquisition functions on banana function.

if ~exist('stats','var') || ~exist('vp','var'); vp = []; stats = []; end

close all;
[vp,stats] = vbmc2020_demo2d(vp,stats);

mypath = '.';
figname = 'vbmc2020_fig1';
% saveas(gcf,[mypath filesep() figname '.pdf']);
% set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;
print([mypath filesep() figname '.pdf'],'-dpdf','-bestfit');
