%VBMC2020_figSpost Plot marginal posteriors

close all;

figname = {'vbmc2020_Spost1','vbmc2020_Spost2'};
probset = 'vbmc20';
algo = 'vbmc@renewdefvarimiqrpluswup5fast';
% algo = 'parallelgp@v3';

for iFig = 1:numel(figname)
    figure(iFig);    
    
    problem = [];
    switch iFig
        case 1
            problem{1} = {'wood2010',1,'D1',{'log(r)','\phi','\sigma_\epsilon'},'Ricker'};
            problem{2} = {'krajbich2010',1,'S1',{'\sigma_\epsilon','d','\beta','\lambda'},'aDDM (S1)'};
            problem{3} = {'krajbich2010',2,'S2',{'\sigma_\epsilon','d','\beta','\lambda'},'aDDM (S2)'};
            problem{4} = {'acerbi2012',1,'S1',{'w_s','w_m','\mu_p','\sigma_p','\lambda'},'Timing'};
            problem{5} = {'acerbidokka2018',1,'S1',{'\sigma_{vis}(c_{low})','\sigma_{vis}(c_{med})','\sigma_{vis}(c_{high})','\sigma_{vest}','\lambda','\kappa'},'Multisensory (S1)'};
        case 2
            problem{1} = {'acerbidokka2018',2,'S2',{'\sigma_{vis}(c_{low})','\sigma_{vis}(c_{med})','\sigma_{vis}(c_{high})','\sigma_{vest}','\lambda','\kappa'},'Multisensory (S2)'};
            problem{2} = {'goris2015b',8,'S8@menoise',{'\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7'},'Neuronal (V1)'};
            problem{3} = {'goris2015b',7,'S7@menoise',{'\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','\theta_7'},'Neuronal (V2)'};
            problem{4} = {'akrami2018b',1,'S1',{'w_L^{(0)}','w_R^{(0)}','w_0','w_L^{(-1)}','w_R^{(-1)}','w_L^{(-2)}','w_R^{(-2)}','w_c','w_s'},'Rodent',[3,9,0.3,-0.4]};
            problem{5} = [];
    end

    nrows = numel(problem);
    ncols = 9;

    N = 5;
    nkde = 2^13;
    fontsize = 24;
    axesfontsize = 18;
    smallfontsize = 14;

    gutter = [.02, .11]; % H V
    margins = [.01 .01 .11 .015]; % L R B T
    hpanel = plotify(nrows,ncols,struct('gutter',gutter,'margins',margins));

    basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';

    for iProb = 1:numel(problem)        
        if isempty(problem{iProb})
            for d = 1:ncols
                axes(hpanel(((iProb-1)*ncols) + d));
                axis off;
            end
            continue;
        end
        prob = problem{iProb}{1};
        id = problem{iProb}{2};

        clear functions;
        probstruct = infprob_init(probset,prob,id,[],1,struct('SkipEval',true));
        LB = probstruct.LB;
        UB = probstruct.UB;

        temp = load([prob '_marginals.mat']);

        bounds = temp.MarginalBounds{id};
        pdf = temp.MarginalPdf{id};
        Nx = size(pdf,2);

        D = size(bounds,2);

        % Load approximate posteriors
        cd([basefolder filesep probset '@' prob filesep problem{iProb}{3}]);

        for ii = 1:N
            try
                temp = load([algo '@' num2str(ii) '.mat']);
                X{ii} = temp.history{1}.Output.post.samples;
            catch
                X{ii} = [];
            end
        end

        for d = 1:D
            axes(hpanel(((iProb-1)*ncols) + d));
            xx = linspace(bounds(1,d),bounds(2,d),Nx);
            plot([xx(1),xx(end)],[0 0],'w-','LineWidth',3); hold on;        
            plot(xx,pdf(d,:),'r-','LineWidth',4);
            hold on;
            box off;
            xmin = bounds(1,d); xmax = bounds(2,d);

            for ii = 1:numel(X)
                [~,yy1,xmesh] = kde1d(X{ii}(:,d),nkde,LB(d),UB(d));
                plot(xmesh,yy1,'-','LineWidth',1,'Color',0*[1,1,1]);
                xmin = min(xmin,xmesh(1));
                xmax = max(xmax,xmesh(end));
            end

            xlim([xmin,xmax]);
            set(gca,'YTick',0,'FontSize',smallfontsize);
            xlabel(problem{iProb}{4}{d},'FontSize',axesfontsize);
        end

        for d = D+1:ncols
            axes(hpanel(((iProb-1)*ncols) + d));
            axis off;        
        end

        if numel(problem{iProb}) > 5
            tt = problem{iProb}{6};
        else
            tt = [iProb,D+1,-0.1,1];
        end

        name = problem{iProb}{5};
        axes(hpanel(((tt(1)-1)*ncols) + tt(2)));
        text(tt(3),tt(4),name,'Units','normalized','FontSize',fontsize,'FontWeight','bold');    

        drawnow;
    end

    % Save figure
    %pos = [1 41 1920 963];
    pos = [1 41 1920 723];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    print([basefolder filesep() figname{iFig} '.pdf'],'-dpdf','-bestfit');    
end

cd(basefolder);