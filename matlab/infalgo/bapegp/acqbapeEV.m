function acq = acqbapeEV(Xs,vbmodel,gp,options,transpose_flag)
%ACQBAPEEV Exponentiated variance acquisition function for BAPE.

if nargin < 5 || isempty(transpose_flag); transpose_flag = false; end

% Threshold on GP variance, try not to go below this
TolVar = options.TolGPVar;

% Transposed input (useful for CMAES)
if transpose_flag; Xs = Xs'; end

% Probability density of vbgmm at test points
% p = vbgmmpdf(vbmodel,Xs')';

% GP mean and variance for each hyperparameter sample
[~,~,fmu,fs2] = gplite_pred(gp,Xs,[],1);

Ns = size(fmu,2);
fbar = sum(fmu,2)/Ns;   % Mean across samples
vbar = sum(fs2,2)/Ns;   % Average variance across samples
if Ns > 1; vf = sum((fmu - fbar).^2,2)/(Ns-1); else; vf = 0; end  % Sample variance
vtot = vf + vbar;       % Total variance

acq = -(2*fbar + vtot + log(expm1(vtot)));

% Regularization: penalize points where GP uncertainty is below threshold
if TolVar > 0
    idx = vtot < TolVar;
    if any(idx)
        acq(idx) = acq(idx) .* exp(-(TolVar./vtot(idx)-1));
    end
end

acq = max(acq,-realmax);

% Transposed output
if transpose_flag; acq = acq'; end

end