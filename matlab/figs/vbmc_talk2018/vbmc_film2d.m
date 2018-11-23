function stats = vbmc_film2d(fun,stats,plotbnd)
%VBMC_DEMO2D Demo movie of VBMC at work (only for 2D problems).

if nargin < 1 || isempty(fun); fun = @rosenbrock_test; end
if nargin < 2 || isempty(stats)
    rng(0);
    [~,~,~,~,~,~,stats] = vbmc(fun,[-1 -1],-Inf,Inf,-3,3);
end
if nargin < 3 || isempty(plotbnd)
    vp = stats.vp(end);
    xrnd = vbmc_rnd(vp,1e6);
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
Npanels = 3;

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

grid = [2 1 3];

Niters = numel(stats.iter);

for iIter = 0:Niters
    
    close all;
    h = plotify(grid,'gutter',[0.075 0.075],'margins',[.05 .02 .15 .075]);

    for iPlot = 1:Npanels
        axes(h(iPlot));
        % subplot(1,3,iPlot);

%         if iIter == 0 && iPlot ~= 2         
%             axis off;
%             continue;
%         end
        
        %[X1,X2] = meshgrid(x1,x2);
        %tmp = cat(2,X2',X1');
        %xx = reshape(tmp,[],2);
        xx = combvec(x1,x2)';

        if iPlot == 1; vpflag = true; else; vpflag = false; end

        elboflag = false;
        if vpflag
            vp = stats.vp(max(1,iIter));
            yy = vbmc_pdf(vp,xx);
            titlestr = ['Iteration ' num2str(iIter)];
            if iIter == 0
                titlestr = [titlestr ' (initial design)'];
            elseif iIter > 1 && ~stats.warmup(iIter) && stats.warmup(iIter-1)
                titlestr = [titlestr ' (end of warm-up)'];
            elseif stats.warmup(iIter)
                titlestr = [titlestr ' (warm-up)'];
            end
        elseif iPlot == 2
            lnyy = zeros(size(xx,1),1);
            for ii = 1:size(xx,1)
                lnyy(ii) = fun(xx(ii,:));
            end
            yy = exp(lnyy);
            Z = sum(yy(:))*dx1*dx2;
            yy = yy/Z;
            titlestr = ['Target density'];
        else
            elboflag = true;
        end

        if elboflag
            idx = 1:iIter;
            iter = stats.iter(idx);
            elbo = stats.elbo(idx);
            elbosd = stats.elboSD(idx);
            beta = 1.96;
            if iIter == 1
                patch([[1 1.2],[1.2 1]],[(elbo + beta*elbosd)*[1 1], (elbo - beta*elbosd)*[1 1]],[1 0.8 0.8],'LineStyle','none'); hold on;                
            elseif iIter > 1
                patch([iter,fliplr(iter)],[elbo + beta*elbosd, fliplr(elbo - beta*elbosd)],[1 0.8 0.8],'LineStyle','none'); hold on;
            end
            if iIter == 0
                hl(1) = plot(1,0,'r','LineWidth',2); hold on;                
            else
                hl(1) = plot(iter,elbo,'r','LineWidth',2); hold on;
            end
            hl(2) = plot([stats.iter(1),stats.iter(end)],log(Z)*[1 1],'k','LineWidth',2);
            if iIter >= 1
                scatter(iter(end),elbo(end),50,'r','FaceColor','r');
            end
            titlestr = 'Model evidence';
            xlim([0.9, stats.iter(end)+0.1]);
            ylims = [floor(min(stats.elbo)-0.5),ceil(max(stats.elbo)+0.5)];
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
            hll = legend(hl,'ELBO','LML');        
            set(hll,'Location',loc,'Box','off');
            box off;
        else
            hold on;

            if vpflag

                % Plot data
                X_train = stats.gp(max(1,iIter)).X;
                    
                if iIter <= 1
                    idx_new = true(size(X_train,1),1);
                else
                    X_trainold = stats.gp(iIter-1).X;
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
            
            if iIter >= 1 || iPlot == 2
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

            xlim([LB(1),UB(1)]);
            ylim([LB(2),UB(2)]);
            set(gca,'TickLength',get(gca,'TickLength')*2);
        end    

        title(titlestr);
        set(gca,'TickDir','out');
    end

    set(gcf,'Color','w');

    pos = [20,20,750,300];
    set(gcf,'Position',pos);
    set(gcf,'Units','inches'); pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    drawnow;
    
    mypath = '.';
    figname = ['demo_' num2str(iIter)];
    saveas(gcf,[mypath filesep() figname '.pdf']);
    
    
%    pause
end
    
end