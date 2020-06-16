function stats = vbmc2020_film2d(fun,stats,plotbnd,options)
%VBMC_DEMO2D Demo movie of VBMC at work (only for 2D problems).

if nargin < 4; options = []; end
if isempty(options)
    options{1}.MaxFunEvals = 150;
    options{1}.MinFunEvals = 150;
    options{1}.MinFinalComponents = 0;
    options{1}.WarpRotoScaling = false;
    options{1}.SearchAcqFcn = @acqfsn2_vbmc;
%    options{1}.SearchAcqFcn = @acqviqr_vbmc;
    options{1}.Text = 'OLD';
    
    options{2} = options{1};
%    options{2}.WarpRotoScaling = true;    
    options{2}.SearchAcqFcn = @acqviqr_vbmc;    
    options{2}.Text = 'NEW';
end

if isstruct(options); options = {options}; end
if isstruct(stats); stats = {stats}; end

Nalgos = max(numel(stats),numel(options));

if nargin < 1 || isempty(fun)
    llfun = @(x) rosenbrock_test(x,2);
    prior_mu = zeros(1,2);
    prior_var = 3^2*ones(1,2);
    lpriorfun = @(x) ...
        -0.5*sum((x-prior_mu).^2./prior_var,2) ...
        -0.5*log(prod(2*pi*prior_var));    
    fun{1} = @(x) lpostfun(x,llfun,lpriorfun);
    fun{2} = @(x) rosenbrock_test(x,0) + lpriorfun(x);
end

if iscell(fun)
    fun_noiseless = fun{2};
    temp = fun{1};
    fun = temp;
else
    fun_noiseless = fun;
end

% Noise test
x0 = [0 0];
y1 = fun(x0); y2 = fun(x0);
noisy_flag = (y1 ~= y2);

if nargin < 2 || isempty(stats)
    for iAlgo = 1:Nalgos
        rng(iAlgo);
        options{iAlgo}.SpecifyTargetNoise = noisy_flag;
        [~,~,~,~,~,~,stats{iAlgo}] = vbmc(fun,[-1 -1],-Inf,Inf,-3,3,options{iAlgo});
    end
end
if nargin < 3 || isempty(plotbnd)
    xrnd = [];
    for iAlgo = 1:Nalgos
        vp = stats{iAlgo}.vp(end);
        xrnd = [xrnd; vbmc_rnd(vp,ceil(1e6/Nalgos))];
    end
    for i = 1:size(xrnd,2)
        LB(i) = floor(quantile1(xrnd(:,i),0.01) - 0.5);
        UB(i) = ceil(quantile1(xrnd(:,i),0.99)+0.5);
    end
else
    LB = plotbnd(1,:);
    UB = plotbnd(2,:);
end

tolx = 1e-3;
Nx = 128;

x1 = linspace(LB(1)+tolx,UB(1)-tolx,Nx);
x2 = linspace(LB(2)+tolx,UB(2)-tolx,Nx);
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);

% Np = 5;
% grid = [];
% for i = 1:(Npanels-2)/2
%     grid = [grid, [i*ones(1,Np); (i+(Npanels-2)/2)*ones(1,Np)]];    
% end
% grid = [grid, [0,Npanels*ones(1,Np);0,(Npanels-1)*ones(1,Np)]];

%labels{1} = 'A';
%labels{Npanels-1} = 'C';
%labels{Npanels} = 'B';

grid = [Nalgos+1, 1:Nalgos, Nalgos+2];
Npanels = numel(grid);

for iAlgo = 1:Nalgos; Niters(iAlgo) = numel(stats{iAlgo}.iter); end
Niters = min(Niters);

for iIter = 0:Niters
    
    close all;
    h = plotify(grid,'gutter',[0.075 0.075],'margins',[.05 .02 .15 .075]);

    for iPlot = 1:Npanels
        axes(h(iPlot));
        
        xx = combvec(x1,x2)';

        % Panels 1-Nalgos are VBMC
        if iPlot <= Nalgos; vpflag = true; else; vpflag = false; end

        elboflag = false;
        if vpflag
            vp = stats{iPlot}.vp(max(1,iIter));
            yy = vbmc_pdf(vp,xx);
            titlestr = ['Iteration ' num2str(iIter)];
            if iIter == 0
                titlestr = [titlestr ' (initial design)'];
            elseif iIter > 1 && ~stats{iPlot}.warmup(iIter) && stats{iPlot}.warmup(iIter-1)
                titlestr = [titlestr ' (end of warm-up)'];
            elseif stats{iPlot}.warmup(iIter)
                titlestr = [titlestr ' (warm-up)'];
            end
            if isfield(options{min(end,iPlot)},'Text') && ~isempty(options{min(end,iPlot)}.Text)
                titlestr = ['VBMC (' options{min(end,iPlot)}.Text ')'];
            end            
        elseif iPlot == Nalgos + 1
            lnyy = zeros(size(xx,1),1);
            for ii = 1:size(xx,1)
                lnyy(ii) = fun_noiseless(xx(ii,:));
            end
            yy = exp(lnyy);
            Z = sum(yy(:))*dx1*dx2;
            yy = yy/Z;
            if noisy_flag
                titlestr = 'Noisy target density';                
            else
                titlestr = 'Target density';
            end
            
            if noisy_flag
                s = contour(x1,x2,reshape(yy',[Nx,Nx])'); hold on;
                lnyy = zeros(size(xx,1),1);
                for ii = 1:size(xx,1)
                    lnyy(ii) = fun(xx(ii,:));
                end
                yy = exp(lnyy);
                Z2 = sum(yy(:))*dx1*dx2;
                yy = yy/Z2;
            end
            
            
        else
            elboflag = true;
        end

        if elboflag
            idx = 1:iIter;
            col = [0 0 1; 1 0 0];
            
            for iAlgo = 1:Nalgos            
                iter = stats{iAlgo}.iter(idx);
                elbo = stats{iAlgo}.elbo(idx);
                elbosd = stats{iAlgo}.elbo_sd(idx);
                beta = 1.96;
                if iIter == 1
%                    patch([[1 1.2],[1.2 1]],[(elbo + beta*elbosd)*[1 1], (elbo - beta*elbosd)*[1 1]],[1 0.8 0.8],'LineStyle','none'); hold on;                
                    patch([[1 1.2],[1.2 1]],[(elbo + beta*elbosd)*[1 1], (elbo - beta*elbosd)*[1 1]],col(iAlgo,:),'FaceAlpha',0.2,'LineStyle','none'); hold on;
                elseif iIter > 1
%                    patch([iter,fliplr(iter)],[elbo + beta*elbosd, fliplr(elbo - beta*elbosd)],[1 0.8 0.8],'LineStyle','none'); hold on;
                    patch([iter,fliplr(iter)],[elbo + beta*elbosd, fliplr(elbo - beta*elbosd)],col(iAlgo,:),'FaceAlpha',0.2,'LineStyle','none'); hold on;
                end
                if iIter == 0
                    hl(iAlgo) = plot(1,0,'-','Color',col(iAlgo,:),'LineWidth',2); hold on;                
                else
                    hl(iAlgo) = plot(iter,elbo,'-','Color',col(iAlgo,:),'LineWidth',2); hold on;
                end
            end
            
            % Plot true log marginal likelihood
            hl(Nalgos+1) = plot([stats{iAlgo}.iter(1),stats{iAlgo}.iter(end)],log(Z)*[1 1],'k','LineWidth',2);
            
            elbo_min = min(stats{iAlgo}.elbo);
            elbo_max = max(stats{iAlgo}.elbo);
            if iIter >= 1
                for iAlgo = 1:Nalgos
                    iter = stats{iAlgo}.iter(idx);
                    elbo = stats{iAlgo}.elbo(idx);
                    elbo_min = min(elbo_min,min(stats{iAlgo}.elbo));
                    elbo_max = max(elbo_max,max(stats{iAlgo}.elbo));
                    scatter(iter(end),elbo(end),50,col(iAlgo,:),'MarkerFaceColor',col(iAlgo,:));
                end
            end
            
            titlestr = 'Model evidence';
            xlim([0.9, stats{iAlgo}.iter(end)+0.1]);
            ylims = [floor(elbo_min-0.5),ceil(elbo_max+0.5)];
            ylim(ylims);
            xxt = sort(unique([1 5:5:iIter iIter]));
            xticks(xxt(xxt <= iIter));
            yticks([ylims(1),round(log(Z),2),ylims(2)])
            xlabel('Iterations');        
            if log(Z) < mean(ylims)
                loc = 'NorthEast';
            else
                loc = 'SouthEast';
            end
            
            for iAlgo = 1:Nalgos
                if isfield(options{min(end,iAlgo)},'Text') && ~isempty(options{min(end,iAlgo)}.Text)
                    legtext{iAlgo} = options{min(end,iAlgo)}.Text;
                else
                    legtext{iAlgo} = 'ELBO';
                end
            end
            legtext{Nalgos+1} = 'LML';
            hll = legend(hl,legtext{:});        
            set(hll,'Location',loc,'Box','off');
            box off;
        else
            hold on;

            if vpflag

                % Plot data
                X_train = stats{iPlot}.gp(max(1,iIter)).X;
                    
                if iIter <= 1
                    idx_new = true(size(X_train,1),1);
                else
                    X_trainold = stats{iPlot}.gp(iIter-1).X;
                    idx_new = false(size(X_train,1),1);
                    [~,idx_diff] = setdiff(X_train,X_trainold,'rows');
                    idx_new(idx_diff) = true;
                end
                idx_old = ~idx_new;

                X_train = warpvars(X_train,'inv',vp.trinfo);
                
                siz = 12;
                if any(idx_old)
                    scatter(X_train(idx_old,1),X_train(idx_old,2),siz,'ok');                            
                end
                if any(idx_new)
                    col = [0.6 0.6 1];
                    scatter(X_train(idx_new,1),X_train(idx_new,2),siz,'o','MarkerFaceColor',col,'MarkerEdgeColor',col);
                end
            end
            
            if iIter >= 1 || iPlot == Nalgos + 1
                s = contour(x1,x2,reshape(yy',[Nx,Nx])');
            end
            
            if vpflag && iIter >= 1
                % Plot component centers
                mu = warpvars(vp.mu','inv',vp.trinfo);
                plot(mu(:,1),mu(:,2),'xr','LineStyle','none');                
            end

            % s.EdgeColor = 'None';
            view([0 90]);
            xlabel('x_1');
            ylabel('x_2');
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            box off;

            xlim([LB(1),UB(1)]);
            ylim([LB(2),UB(2)]);
            set(gca,'TickLength',get(gca,'TickLength')*2);
        end    

        title(titlestr);
        set(gca,'TickDir','out');
    end

    set(gcf,'Color','w');

    pos = [20,20,750/3*(Nalgos+2),300];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    
    if 0
        mypath = '.';
        figname = ['demo_' num2str(iIter)];
        saveas(gcf,[mypath filesep() figname '.pdf']);
    else
        endIters = max(3,floor(Niters/2));
        f = getframe(gcf);
        if iIter == 0
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,Niters+endIters+1) = 0;
        else
            im(:,:,1,iIter+1) = rgb2ind(f.cdata,map,'nodither');
        end
        if iIter == Niters
            if Nalgos == 1
                axes(h(1));    
                titlestr = ['Iteration ' num2str(iIter) ' (stable)'];
                title(titlestr);
                drawnow;
            end
            f = getframe(gcf);
            for iFinal = 1:endIters
                im(:,:,1,Niters+1+iFinal) = rgb2ind(f.cdata,map,'nodither');
            end
            filename = 'vbmc2020_demo';
            imwrite(im,map,[filename '.gif'],'DelayTime',0,'LoopCount',inf);
        end
    end
        
    
    
%    pause
end
    
end