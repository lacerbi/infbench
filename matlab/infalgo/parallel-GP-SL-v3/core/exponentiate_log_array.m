function [res,maxu] = exponentiate_log_array(log_vals)
% Computes exp(log_vals) up to multiplicative normalisation factor i.e. output is scaled
% so that the output does not (completely) under- or overflow.

if any(isnan(log_vals))
    warning('NaN values encountered');
end

maxu = max(log_vals);
if exp(maxu) == Inf || exp(maxu) == 0
    % rescale so that maximum is exp(0)==1
    res = exp(log_vals - maxu);
else
    % scaling was not necessary
    res = exp(log_vals);
    maxu = 0;
end
end


