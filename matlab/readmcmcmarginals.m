function [MarginalBounds,MarginalPdf] = readmcmcmarginals(prob,subprobs)
%READMCMCMARGINALS Read marginals from MCMC runs.

% Example: [MarginalBounds,MarginalPdf] = readmcmcmarginals('goris2015',7:8)
% (To be launched from the goris2015mcmc data folder).

nkde = 2^13;

for iSubprob = 1:numel(subprobs)
    subprob = subprobs(iSubprob);    
    X = mergesamples([prob '_mcmc_n' num2str(subprob) '*'],Inf,0,0);
    D = size(X,2);
    
    MarginalBounds{subprob} = zeros(2,D);
    MarginalPdf{subprob} = zeros(D,nkde);
    
    for i = 1:D
        [~,yy1,xmesh] = kde(X(:,i),nkde);
        MarginalBounds{subprob}(:,i) = [xmesh(1);xmesh(end)];
        yy1 = yy1/(qtrapz(yy1)*(xmesh(2)-xmesh(1))); % Ensure normalization
        MarginalPdf{subprob}(i,:) = yy1;
    end    
end

savefolder = fileparts(which(['infbench_' prob '.m']));
save([savefolder filesep prob '_marginals.mat'],'MarginalBounds','MarginalPdf');
