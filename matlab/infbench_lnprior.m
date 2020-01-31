function lnp = infbench_lnprior(X,probstruct)
%INFBENCH_LNPRIOR Compute standard log prior for inference benchmark.

if ~isfield(probstruct,'PriorType') || isempty(probstruct.PriorType)
    PriorType = 'gaussian';
else
    PriorType = probstruct.PriorType;
end

switch lower(PriorType)    
    case 'gaussian'
        lnp = -0.5*sum(log(2*pi*probstruct.PriorVar),2) ...
            -0.5*sum(bsxfun(@rdivide, ...
            bsxfun(@minus,X,probstruct.PriorMean).^2, ...
            probstruct.PriorVar),2);        
    case 'uniform'        
        lnp = -log(probstruct.PriorVolume)*ones(size(X,1),1);        
    otherwise
        error('Unknown prior for INFBENCH.');
end