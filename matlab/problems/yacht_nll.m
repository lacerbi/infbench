function y = yacht_nll(x,infprob,priorflag)

if nargin < 3 || isempty(priorflag); priorflag = 1; end

y = -infbench_yacht(x,infprob);

if priorflag
    infprob.PriorMean = infprob.Prior.Mean;
    infprob.PriorVar = diag(infprob.Prior.Cov)';
    lnp = infbench_lnprior(x(:)',infprob);
    y = y - lnp;
end
    
end