function output = infbench_results(prob_list,algo_list,basefolder)
%INFBENCH_PERFORMANCE Collect benchmark results.

if nargin < 3; basefolder = []; end

if isempty(basefolder)
    % Default base folder for benchmark output files
%    basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
    basefolder = 'C:\Users\luigi\Dropbox\scratch\infbench';
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
                
        output.lnZerr{iProb,iAlgo} = [];
        output.gsKL{iProb,iAlgo} = [];
        output.MMTV{iProb,iAlgo} = [];
        
        for iRun = 1:Nmax
           filename =  [algo_list{iAlgo} '@' num2str(iRun) '.mat'];
           try
               temp = load(filename);
           catch
               temp = [];
           end
           if isempty(temp); continue; end

           for id = 1:numel(temp.history)
               history = temp.history{id};
               
               if isfield(temp,'speedtest')
                   t = temp.speedtest.start(:,1:4) + temp.speedtest.end(:,1:4);
                   history.speedtest = sum(t(:));
               else
                   history.speedtest = NaN;
               end
               
               % Write operations here
               lnZerr = abs(history.Output.post.lnZ - history.lnZpost_true);
               gsKL = history.Output.post.gsKL;
               MMTV = mean(history.Output.post.gsKL);

               output.lnZerr{iProb,iAlgo} = [output.lnZerr{iProb,iAlgo}; lnZerr];
               output.gsKL{iProb,iAlgo} = [output.gsKL{iProb,iAlgo}; gsKL];
               output.MMTV{iProb,iAlgo} = [output.MMTV{iProb,iAlgo}; MMTV];               
           end           
        end
    end
end