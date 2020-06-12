%VBMC2020_figSpostcorner Cornerplot of posteriors

close all;
save_flag = true;

figname = {'vbmc2020_Spostcorner1','vbmc2020_Spostcorner2'};
probset = 'vbmc20';
algo = 'vbmc@renewdefvarimiqrpluswup5fast';
idx = 101;
txtpanel = {'A','B'};
txtpanel2 = {'True','VBMC-VIQR'};

fontsize = 24;
axesfontsize = 18;
smallfontsize = 14;

problem = {'acerbi2012',1,'S1',{'w_s','w_m','\mu_p','\sigma_p','\lambda'},'Timing'};

basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';

prob = problem{1};
id = problem{2};

clear functions;
probstruct = infprob_init(probset,prob,id,[],1,struct('SkipEval',true));
LB = probstruct.LB;
UB = probstruct.UB;

% Load ground-truth posterior
X = [];
cd([basefolder filesep prob 'mcmc']);
X{1} = mergesamples([prob '_mcmc_n' num2str(id) '_*'],[],0,0);

D = size(X{1},2);

% Load approximate posteriors
cd([basefolder filesep probset '@' prob filesep problem{3}]);

temp = load([algo '@' num2str(idx) '.mat']);
X{2} = temp.history{1}.Output.post.samples;

LB = min([X{1};X{2}]);
UB = max([X{1};X{2}]);

for iFig = 1:numel(figname)

%    gutter = [.02, .11]; % H V
%    margins = [.01 .01 .11 .015]; % L R B T
%    hpanel = plotify(nrows,ncols,struct('gutter',gutter,'margins',margins));

    %margins = [.1 .01 .12 .01];
    margins = [.1 .01 .12 .05]; % L R B T

    [~,ax] = cornerplot(X{iFig},problem{4},[],[LB;UB],0,margins);
    for iAx = 1:numel(ax)
        if isnan(ax(iAx)); continue; end        
        axes(ax(iAx));
        xl = get(gca,'Xlabel'); xtxt = xl.String;
        yl = get(gca,'Ylabel'); ytxt = yl.String;        
        set(gca,'FontSize',smallfontsize);
        xlabel(xtxt,'FontSize',axesfontsize);
        ylabel(ytxt,'FontSize',axesfontsize);        
    end
    
    axes(ax(1));    
    text(-0.55,1.2,txtpanel{iFig},'Units','normalized','FontSize',fontsize,'FontWeight','bold');    

    text(3.5,0.5,txtpanel2{iFig},'Units','normalized','FontSize',axesfontsize);    
    
    drawnow;

    % Save figure
        %pos = [1 41 1920 963];
    pos = [73   210   871   674];
    set(gcf,'Position',pos);
    if save_flag
        set(gcf,'Units','inches'); pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Manual','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        drawnow;
        print([basefolder filesep() figname{iFig} '.pdf'],'-dpdf','-bestfit');
    end
end

cd(basefolder);