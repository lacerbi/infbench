function vbmc20_performance_jitternoise(prob,subprob,algo_list,metric_idx,Nmax)
%VBMC20_PERFORMANCE_JITTERNOISE Performance with added jitter to estimated noise

if nargin < 1; prob = []; end
if nargin < 2; subprob = []; end
if nargin < 3; algo_list = []; end
if nargin < 4 || isempty(metric_idx); metric_idx = 2; end
if nargin < 5 || isempty(Nmax); Nmax = 50; end  % Max number of output files

probset = 'vbmc20';
if isempty(prob); prob = 'acerbi2012'; end
if isempty(subprob); subprob = 'S101'; end

if isempty(algo_list)
    algos_list{1} = {'vbmc@renewdefvarimiqrpluswup5fast'};
    % algos_list{1} = {'vbmc@renewdefimiqrpluswup5noacq','vbmc@renewdefmipluswup4gpsvp','vbmc@renewdefimiqrplus5longvpgps','vbmc@renewdefvarimiqrpluswup5fast'};
    algos_list{2} = {'vbmc@renewdefvarimiqrpluswup5fast'};
    algos_list{3} = {'vbmc@renewdefvarimiqrpluswup5fast'};
    algo_list = algos_list{metric_idx};
end

fontsize = 18;
axesfontsize = 14;
bootstrap = 1e4;

thresh = [1,0.2,1];

noise = 2;
%jitter_levels = [0,0.2,0.4,0.5,0.7];
jitter_levels = 0:0.1:0.7;
for iNoise = 1:numel(jitter_levels)
    if jitter_levels(iNoise) == 0
        sub = [subprob '@n' num2str(noise) 'noise'];
    else
        sub = [subprob '@n' num2str(noise) 'j' num2str(jitter_levels(iNoise)) 'noise'];
    end
    prob_list{iNoise} = {probset,prob,sub};
end

basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
out = infbench_collect(prob_list,algo_list,basefolder,Nmax);

ytick = [0.0001,0.001,0.01,0.1,1,10,1e2,1e3,1e4];
yticklabel = {'10^{-4}','10^{-3}','0.01','0.1','1','10','100','10^3','10^4'};

switch metric_idx
    case 1; ylims = [0.01,1e3]; ystring = 'LML loss'; logy_flag = true;
    case 2; ylims = [0,1]; ystring = 'MMTV'; ytick = 0:0.2:1; yticklabel = []; logy_flag = false;
    case 3; ylims = [0.01,1e3]; ystring = 'gsKL'; logy_flag = true;
end

for iAlgo = 1:numel(algo_list)
    
    defaults = infbench_defaults('style',[],[],algo_list{iAlgo});
    algo_name{iAlgo} = defaults.name;
    
    xx = jitter_levels;
    for iNoise = 1:numel(jitter_levels)
        yy = out{iNoise,iAlgo}(:,metric_idx);
        y(iNoise) = median(yy);        
        if bootstrap > 0
            yy_boot = bootstrp(bootstrap,@median,yy);
            y_up(iNoise) = quantile(yy_boot,0.975,1);
            y_down(iNoise) = quantile(yy_boot,0.025,1);
        else
            y_up(iNoise) = quantile(yy,0.75);
            y_down(iNoise) = quantile(yy,0.25);
        end
    end
    style = [defaults.linestyle];
    lincol = defaults.color;
    lw = defaults.linewidth;
    xxerr = [xx, fliplr(xx)];
    yyerr = [y_down, fliplr(y_up)];
    fill(xxerr,yyerr,lincol,'FaceAlpha',0.5,'LineStyle','none'); hold on;
    h(iAlgo) = plot(xx,y,style,'Color',lincol,'LineWidth',lw);
    
    xlims = [jitter_levels(1),jitter_levels(end)];
    plot(xlims,thresh(metric_idx)*[1 1],'k--','Linewidth',0.5);
    
end

hl = legend(h,algo_name{:});
set(hl,'Box','off','Location','NorthWest','FontSize',fontsize);


set(gca,'Xlim',xlims,'Ylim',ylims,'YTick',ytick);
if logy_flag
    set(gca,'Yscale','log','YMinorTick','off');
end
if ~isempty(yticklabel); set(gca,'YTickLabel',yticklabel); end
set(gca,'TickDir','out','TickLength',3*get(gca,'TickLength')); 
set(gca,'XTick',jitter_levels,'FontSize',axesfontsize);
box off;
xlabel('Log-likelihood noise \sigma_{obs}','FontSize',fontsize);
ylabel(ystring,'FontSize',fontsize);
set(gcf,'Color','w');
drawnow;

end


%--------------------------------------------------------------------------
function output = infbench_collect(prob_list,algo_list,basefolder,Nmax)
%INFBENCH_COLLECT Template function for collecting benchmark outputs.

if nargin < 3; basefolder = []; end

if isempty(basefolder)
    % Default base folder for benchmark output files
    basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
end

BaseSpeedTest = 13.7660; % Laptop speed

clear S;
for iProb = 1:numel(prob_list)
    
    probset = prob_list{iProb}{1};
    prob = prob_list{iProb}{2};
    subprob = prob_list{iProb}{3};
    
    folder = [basefolder filesep probset '@' prob filesep subprob];
    cd(folder);
    
    for iAlgo = 1:numel(algo_list)
        
        probstr = [probset '@' prob '@' subprob];
        fprintf('Problem: %s. Algorithm: %s...\n',probstr,algo_list{iAlgo});
                
        output{iProb,iAlgo} = [];
        
        for iRun = 1:Nmax
           filename =  [algo_list{iAlgo} '@' num2str(iRun) '.mat'];
           try
               temp = load(filename);
               history = temp.history{1};               
                if isfield(temp,'speedtest')
                    t = temp.speedtest.start(:,1:4) + temp.speedtest.end(:,1:4);
                    history.speedtest = sum(t(:));
                else
                    history.speedtest = NaN;
                end               
           catch
               temp = [];
           end
           if isempty(temp); continue; end
           
           idx_valid = history.SaveTicks <= history.TotalMaxFunEvals;
           last = find(idx_valid,1,'last');
           speedfactor = BaseSpeedTest/history.speedtest;                       
                      
           % Write operations here
           
           FuncCumTime = cumsum(history.FuncTime);
           overhead_new = speedfactor*(history.ElapsedTime(last) - FuncCumTime(last))./history.SaveTicks(last);
       
           lnZ_err = abs(history.lnZpost_true - history.Output.post.lnZ);
           MMTV = mean(history.Output.post.MTV);
           gsKL = history.Output.post.gsKL;
           
           output_new = [lnZ_err,MMTV,gsKL];
           output{iProb,iAlgo} = [output{iProb,iAlgo}; output_new];
        end
    end
end

end