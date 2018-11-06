function gpplotdemo(seed)

if nargin < 1 || isempty(seed); seed = 1; end

mypath = fileparts(mfilename('fullpath'));

rng(seed);
sn = 1e-3;
hyp{1} = [-0.3;0;log(sn);0];
hyp{2} = [0.7;0;log(sn);0];
hyp{3} = [-0.3;0;log(sn);-2;0;0.5];

meanfun = {'const','const','quad'};
titlestr = {'short length scale';'long length scale';'quadratic mean fcn.'};

figname = {'gpdemo_prior','gpdemo_post'};

fontsize = 18;
axesfontsize = 14;
close all;

for iFig = 1:2
    figure(iFig);
    for iPlot = 1:3
        subplot(1,3,iPlot);
        if iFig == 1
            X = zeros(0,1);
            y = [];
        else
            X = [-2; 0; 2];
            y = [0; -1; 1];
        end
        gp = gplite_post(hyp{iPlot},X,y,meanfun{iPlot});
        Xstar = linspace(-3,3,201)';

        [~,~,fmu,fs2] = gplite_pred(gp,Xstar);
        % plot(xx,fmu,'-k','LineWidth',1); hold on;    
        xx = [Xstar',fliplr(Xstar')];
        fill(xx,[fmu',fliplr(fmu')]+1.96*[sqrt(fs2'),-fliplr(sqrt(fs2'))],0.85*[1 1 1],'LineStyle','none'); hold on;

        for iGP = 1:3
            Ystar = gplite_rnd(gp,Xstar);
            scatter(Xstar,Ystar,'.');
            hold on;
        end
        if size(X,1) > 1
            scatter(X,y,'ko');        
        end
        title(titlestr{iPlot},'FontSize',fontsize);
        xlim([min(Xstar),max(Xstar)]);
        ylim([-5,5]);
        set(gca,'TickDir','out','Xtick',[min(Xstar),0,max(Xstar)],'Ytick',[-4,0,4],'FontSize',axesfontsize);
        xlabel('input, x','FontSize',fontsize);
        if iPlot == 1
            ylabel('output, f(x)','FontSize',fontsize);
        end
        box off;
    end
    set(gcf,'Color','w');

    % Save figure
    pos = [168 519 1400 420];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    saveas(gcf,[mypath filesep() figname{iFig} '.pdf']);
end