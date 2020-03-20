function [nLL,output] = ibl2020_nllfun(theta,params,data)
%NLLFUN Negative log likelihood function

if isempty(theta)
    theta = params.theta;
end

% Assign parameter vector to parameter struct
params1 = ibl2020_setup_params(theta,params);

if nargout == 1
    nLL = params1.model_nLLfun(params1,data);
else
    [nLL,output] = params1.model_nLLfun(params1,data);
end

nLL = sum(nLL);

end