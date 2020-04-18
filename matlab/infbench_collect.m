function output = infbench_collect(prob_list,algo_list,basefolder)
%INFBENCH_COLLECT Template function for collecting benchmark outputs.

if nargin < 3; basefolder = []; end

if isempty(basefolder)
    % Default base folder for benchmark output files
    basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
end

BaseSpeedTest = 13.7660; % Laptop speed
Nmax = 100; % Max number of output files

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
       
           output{iProb,iAlgo} = [output{iProb,iAlgo}; overhead_new];
        end
    end
end