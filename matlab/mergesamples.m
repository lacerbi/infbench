function [X,lls,funccount,Neffs] = mergesamples(filepattern,MaxIter,Nml,psrf_flag)
%MERGESAMPLES Merge MCMC samples from different files.

% Example usage: from the data/acerbidokka2018mcmc folder
% [X,lls,funccount] = mergesamples('acerbidokka2018_mcmc_n1_*',[],2e4,1);

if nargin < 2 || isempty(MaxIter); MaxIter = Inf; end
if nargin < 3 || isempty(Nml); Nml = 2e4; end
if nargin < 4 || isempty(psrf_flag); psrf_flag = true; end
if nargin < 5; infprob = []; end

files = dir(filepattern);
M = numel(files);

Xs = [];
lls = [];
funccount = NaN(M,1);
iter = 1;

for iFile = 1:M
    if files(iFile).isdir; continue; end
    filename = files(iFile).name;
    try
        temp = load(filename);
        
        if isempty(Xs) || isempty(lls)
            N = size(temp.Xs,1);
            D = size(temp.Xs,2);
            Xs = NaN(N,D,M);
            lls = NaN(N,M);
        end
        
        Xs(:,:,iter) = temp.Xs;
        lls(:,iter) = temp.lls(:,1);        
        if isfield(temp.output,'funccount')
            funccount(iter) = temp.output.funccount;
        elseif isfield(temp.output,'nsimu')
            funccount(iter) = temp.output.nsimu;
        else
            funccount(iter) = NaN;
        end
        if isfield(temp.output,'Neff_preburn')
            Neffs.preburn(iter,:) = temp.output.Neff_preburn;
            Neffs.afterburn(iter,:) = temp.output.Neff_afterburn;
        end
        
        fprintf('%d..', iter);
        iter = iter + 1;
        if iter > MaxIter
            fprintf('\nReached maximum number of files %d.\n', MaxIter);
            break; 
        end
    catch
        fprintf('\nCould not read data from file %s.\n', filename);
    end
end
fprintf('\n');

imax = iter-1;

Xs = Xs(:,:,1:imax);
lls = lls(:,1:imax);
funccount = funccount(1:imax);

% Compute PSRF
if psrf_flag
    [R,Neff] = psrf(Xs);
else
    R = NaN; Neff = NaN;
end
fprintf('R_max = %.3f. Ntot = %d. Neff_min = %.1f. Total funccount = %d.\n',max(R),N*M,min(Neff),sum(funccount));

% Reshape X
[N,D,M] = size(Xs);
X = NaN(N*M,D);
for m = 1:M
    X((1:N)+N*(m-1),:) = Xs(:,:,m);
end
fprintf('\n\tMean_mcmc = %s;\n\tCov_mcmc = %s;\n', mat2str(mean(X,1)), mat2str(cov(X)));

if Nml > 0    
    N = size(X,1);
    % Compute approximation of marginal likelihood
    if N > Nml
        idx = round(linspace(1,N,Nml))';
    else
        idx = (1:N)';
    end
    lnZ = vbgmmnormconst('rlr',X(idx,:),lls(idx));
    fprintf('\tlnZ_mcmc = %s;\n', mat2str(lnZ));
end

lls = lls(:);


