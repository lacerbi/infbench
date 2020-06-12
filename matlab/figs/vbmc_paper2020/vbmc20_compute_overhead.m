%VBMC20_COMPUTE_OVERHEAD

prob_list = [];
prob_list{1} = {'vbmc20','wood2010','D1'};
% prob_list{2} = {'vbmc20','price2018','D1'};
prob_list{2} = {'vbmc20','krajbich2010','S1'};
prob_list{3} = {'vbmc20','krajbich2010','S2'};
prob_list{4} = {'vbmc20','acerbi2012','S1'};
prob_list{5} = {'vbmc20','acerbidokka2018','S1'};
prob_list{6} = {'vbmc20','acerbidokka2018','S2'};
prob_list{7} = {'vbmc20','goris2015b','S108@menoise'};
prob_list{8} = {'vbmc20','goris2015b','S107@menoise'};
prob_list{9} = {'vbmc20','akrami2018b','S1'};

algo_list{1} = 'vbmc@renewdefvarimiqrpluswup5fast';
% algo_list{1} = 'vbmc@viqrnogplate';
algo_list{2} = 'vbmc@renewdefmipluswup4gpsvp';
algo_list{3} = 'vbmc@renewdefimiqrplus5longvpgps';
algo_list{4} = 'parallelgp@v3';

algo_name{1} = 'vbmc-viqr';
algo_name{2} = 'vbmc-eig';
algo_name{3} = 'vbmc-imiqr';
algo_name{4} = 'gp-imiqr';

basefolder = 'C:\Users\Luigi\Dropbox\Postdoc\VBMC\data';
out = infbench_collect(prob_list,algo_list,basefolder);

clear X;

for iAlgo = 1:numel(algo_list)
    X{1,iAlgo} = out{1,iAlgo};
    X{2,iAlgo} = [out{2,iAlgo}; out{3,iAlgo}];
    X{3,iAlgo} = out{4,iAlgo};
    X{4,iAlgo} = [out{5,iAlgo}; out{6,iAlgo}];
    X{5,iAlgo} = [out{7,iAlgo}; out{8,iAlgo}];
    X{6,iAlgo} = out{9,iAlgo};
end

for iAlgo = 1:numel(algo_list)
    fprintf('\\texttt{%s} ',algo_name{iAlgo});
    for iProb = 1:size(X,1)
%        fprintf('& %.1f ',mean(X{iProb,iAlgo}));
%        fprintf('& %.1f [%.1f, %.1f]',mean(X{iProb,iAlgo}),quantile(X{iProb,iAlgo},0.025),quantile(X{iProb,iAlgo},0.975));
        fprintf('& $%.1f \\pm %.1f$ ',mean(X{iProb,iAlgo}),std(X{iProb,iAlgo}));
%        fprintf('& %.1f [%.2f]',mean(X{iProb,iAlgo}), ...
%            max(abs(1 - quantile(X{iProb,iAlgo},[0.025,0.975])/mean(X{iProb,iAlgo}))));
    end
    fprintf('\\\\\n');
end
