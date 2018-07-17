function [history,post] = StoreAlgoResults(probstruct,post,Niter,X,y,mu,vvar,Xiter,yiter,TotalTime)
%STOREALGORESULTS Store results of running an inference algorithm.

history = infbench_func(); % Retrieve history
% history.scratch.output = output;
history.TotalTime = TotalTime;
history.Output.X = X;
history.Output.y = y;
if ~probstruct.AddLogPrior      % y stores log posteriors, so add prior now
    lnp = infbench_lnprior(history.Output.X,probstruct);
    history.Output.y = history.Output.y + lnp;
end

if ~isempty(mu)
    post.lnZ = mu(end);
    post.lnZ_var = vvar(end);
    compute_lnZ = false;
else
    compute_lnZ = true;
end
X_train = history.Output.X;
y_train = history.Output.y;
[post.gsKL,post.Mean,post.Cov,lnZ,lnZ_var,post.Mode] = ComputeAlgoStats(X_train,y_train,probstruct,compute_lnZ);
if compute_lnZ
    post.lnZ = lnZ;
    post.lnZ_var = lnZ_var;
end

% Return estimate, SD of the estimate, and gauss-sKL with true moments
if isempty(Niter)
    Niter = find(size(X,1) == history.SaveTicks,1);
end
N = history.SaveTicks(1:Niter);
history.Output.N = N(:)';
if ~isempty(mu)
    history.Output.lnZs = mu(:)';
    history.Output.lnZs_var = vvar(:)';
else
    history.Output.lnZs = NaN(1,Niter);
    history.Output.lnZs_var = NaN(1,Niter);
end

for iIter = 1:Niter
    if isempty(Xiter) || isempty(yiter)
        X_train = history.Output.X(1:N(iIter),:);
        y_train = history.Output.y(1:N(iIter));
    else
        X_train = Xiter{iIter};
        y_train = yiter{iIter};
    end
    [gsKL,Mean,Cov,lnZ,lnZ_var,Mode] = ComputeAlgoStats(X_train,y_train,probstruct,compute_lnZ);
    if compute_lnZ
        history.Output.lnZs(iIter) = lnZ;
        history.Output.lnZs_var(iIter) = lnZ_var;
    end
    history.Output.Mean(iIter,:) = Mean;
    history.Output.Cov(iIter,:,:) = Cov;
    history.Output.gsKL(iIter) = gsKL;
    history.Output.Mode(iIter,:) = Mode;    
end


end