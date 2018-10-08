% FIGURE S1 for VBMC paper. Plot examples of synthetic functions.

close all;
grid = [1 3 5; 2 4 6];
h = plotify(grid,'gutter',[0.1 0.04],'margins',[.05 .02 .075 .05]); % ,'labels',labels

for iPlot = 1:3
    
    rng(0);

    switch iPlot
        case 1
            prob = 'lumpy';
            plotlb = [-0.5,-0.5];
            plotub = [0.5,0.5];
            titlestr = 'Lumpy';
        case 2
            prob = 'studentt';
            plotlb = [-0.5,-0.5];
            plotub = [0.5,0.5];
            titlestr = 'Student';
        case 3
            prob = 'cigar';
            plotlb = [-0.5,-0.5];
            plotub = [0.5,0.5];
            titlestr = 'Cigar';
    end

    clear functions;
    probstruct = infprob_init('vbmc18',prob,2,[],1,[]);
    probstruct.AddLogPrior = true;

    x0 = probstruct.InitPoint;
    PLB = probstruct.PLB;
    PUB = probstruct.PUB;
    LB = probstruct.LB;
    UB = probstruct.UB;
    fun = @(x) infbench_func(x,probstruct);

    if 1
        [vp,~,~,~,output,~,stats] = vbmc(fun,x0,LB,UB,PLB,PUB);
        gp = stats.gp(output.bestiter);
    end
    
    axes(h(iPlot*2-1));
    vbmc_plot2d(fun,plotlb,plotub);
    axis square;
    title(titlestr);
    xlabel('');
    text(-0.4,0.4,'True');
    
    gp = [];
    axes(h(iPlot*2));
    vbmc_plot2d(vp,plotlb,plotub,gp,0);
    axis square;
    text(-0.4,0.4,'VBMC');
    text(0.4,-0.4,['Iteration ' num2str(output.bestiter)],'HorizontalAlignment','right');
    
    drawnow;

end

pos = [20,20,900,500];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;

mypath = '.';
figname = 'vbmc_figS1';
saveas(gcf,[mypath filesep() figname '.pdf']);

