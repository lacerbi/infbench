function [vp,stats] = vbmc2020_demo2d(vp,stats)
%VBMC2020_DEMO2D

if nargin < 1; vp = []; stats = []; end

sigman = 1;
fun = @(x) targetfun(x,sigman);
x0 = [-1 -1];
PLB = [-3 -3];
PUB = [3 3];
plotbnd = [-3,-2; 3,6];

% Set basic VBMC options
options = vbmc('all');
options.ActiveSampleGPUpdate = true;
options.ActiveSampleVPUpdate = true; 
options.WarpRotoScaling = true;
options.UncertaintyHandling = 'on';
options.SpecifyTargetNoise = 1;
options.MaxFunEvals = 100;
options.MinFunEvals = 100;
options.MinFinalComponents = 0;
% options.NSgpMaxMain = 0;

text = {'NPRO','VIQR'};
AcqFuns = {@acqfsn2_vbmc,@acqviqr_vbmc};
% AcqFuns = {@acqmaxiqrpro_vbmc,@acqvarimiqr_vbmc}; options.Plot = 1;
colors = [178,223,138; 0 0 0]/255;

%AcqFuns = {@acqfsn2_vbmc,@acqmi_vbmc,@acqimiqr_vbmc,@acqvarimiqr_vbmc};

if isempty(vp) || isempty(stats)
    for iAcq = 1:numel(AcqFuns)
        rng(iAcq);
        options.SearchAcqFcn = AcqFuns{iAcq};    
        [vp{iAcq},~,~,~,~,~,stats{iAcq}] = vbmc(fun,x0,-Inf,Inf,PLB,PUB,options);
        % xx = vbmc_rnd(vp{iAcq},1e6);
        % cornerplot(xx);
    end
end

LB = plotbnd(1,:);
UB = plotbnd(2,:);

tolx = 1e-3;
Nx = 128;
Npanels = 4;

x1 = linspace(LB(1)+tolx,UB(1)-tolx,Nx);
x2 = linspace(LB(2)+tolx,UB(2)-tolx,Nx);
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);

Np = 5;
grid = [ones(1,Np),0];
for i = 2:(Npanels-1); grid = [grid, i*ones(1,Np)]; end
grid = [grid, [0,Npanels*ones(1,Np)]];

% grid = [reshape(1:Npanels-2,[(Npanels-2)/2,2])',[Npanels;Npanels-1]];
labels{1} = 'A';
labels{2} = 'B';
%labels{3} = 'C';
labels{4} = 'C';

h = plotify(grid,'gutter',[0.05 0.15],'margins',[.05 .02 .175 .1],'labels',labels);

for iPlot = 1:Npanels
    axes(h(iPlot));
    
    xx = combvec(x1,x2)';
        
    vpflag = false;
    elboflag = false;
    if iPlot > 1 && iPlot < Npanels
        vpflag = true;
        vp1 = vp{iPlot-1};
        yy = vbmc_pdf(vp1,xx);
        titlestr = ['VBMC+' text{iPlot-1}];
    elseif iPlot == 1
        lnyy = zeros(size(xx,1),1);
        for ii = 1:size(xx,1)
            lnyy(ii) = targetfun(xx(ii,:),0);
        end
        yy = exp(lnyy);
        Z = sum(yy(:))*dx1*dx2;
        yy = yy/Z;
        titlestr = ['True posterior'];
    else
        elboflag = true;
    end
    
    if elboflag
        
        for iStats = 1:numel(stats)
            stats1 = stats{iStats};
            iter = stats1.funccount;
            elbo = stats1.elbo;
            elbo_sd = stats1.elbo_sd;
            beta = 1.96;
            col = colors(iStats,:);
            colb = 0.8 + col*0.2;
            %patch([iter,fliplr(iter)],[elbo + beta*elbo_sd, fliplr(elbo - beta*elbo_sd)],colb,'LineStyle','none'); hold on;
            patch([iter,fliplr(iter)],[elbo + beta*elbo_sd, fliplr(elbo - beta*elbo_sd)],col,'LineStyle','none','FaceAlpha',0.2); hold on;
            hl(iStats) = plot(iter,elbo,'LineWidth',2,'Color',col); hold on;
        end
        hl = [hl(2),hl(1)];
        % hl(numel(stats)+1) = 
        plot([iter(1),iter(end)],log(Z)*[1 1],'k--','LineWidth',1);
        titlestr = 'Model evidence';
        xlim([9, iter(end)+1]);
        ylims = [-4,-1];
        % ylims = [floor(min(elbo)-0.5),ceil(max(elbo)+0.5)];
        ylim(ylims);
        xticks(20:20:200);
        yticks([ylims(1),round(log(Z),2),ylims(2)])
        xlabel('Likelihood evaluations');        
        if log(Z) < mean(ylims)
            loc = 'NorthEast';
        else
            loc = 'SouthEast';
        end
        hll = legend(hl,text{2},text{1});
        % set(hll,'Location',loc,'Box','off');        
        set(hll,'Units','normalized','Box','off');        
        set(hll,'Position',[0.785 0.2 0.3 0.15]);
    else
        s = contour(x1,x2,reshape(yy',[Nx,Nx])');

        if vpflag
            % Plot component centers
            mu = warpvars_vbmc(vp1.mu','inv',vp1.trinfo);
            hold on;
            plot(mu(:,1),mu(:,2),'xr','LineStyle','none');

            % Plot data
            X = warpvars_vbmc(vp1.gp.X,'inv',vp1.trinfo);
            plot(X(:,1),X(:,2),'.k','LineStyle','none');
        end

        % s.EdgeColor = 'None';
        view([0 90]);
        xlabel('\theta_1');
        ylabel('\theta_2');
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        xlim([LB(1),UB(1)]);
        ylim([LB(2),UB(2)]);
        set(gca,'TickLength',get(gca,'TickLength')*2);
    end    
    
    title(titlestr);
    set(gca,'TickDir','out');
end

set(gcf,'Color','w');

pos = [20,20,900,225];
set(gcf,'Position',pos);
set(gcf,'Units','inches'); pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
drawnow;

end



function [y,y_sd] = targetfun(x,sigman)

sigma2 = 9;     % Prior variance
lprior = -0.5*sum(x.^2,2)/sigma2 - 0.5*size(x,2)*log(2*pi*sigma2);
llike = -sum((x(:,1:end-1) .^2 - x(:,2:end)) .^ 2 + (x(:,1:end-1)-1).^2/100,2);

y = lprior + llike + sigman*randn(size(x,1),1);
y_sd = sigman;

end