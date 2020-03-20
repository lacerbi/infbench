function params = ibl2020_params_new(model_name,data)
%PARAMS_NEW Return new parameter struct for a given model.

if isempty(model_name); model_name = 'exponential'; end

params.model_name = model_name;
params.model_startingfit = [];  % Model used to initialize fits
extra_features = [];

%% Properties of base models

idx = find(model_name=='_',1);
if isempty(idx); idx = numel(model_name)+1; end
base_model = model_name(1:idx-1);
extra_model = model_name(idx:end);

switch base_model
    case 'exponential'
        params.model_nLLfun = @exponential_nll;
        params.model_desc = 'Exponential-averaging observer';
    otherwise
        error('Unknown base model.');
end        

%% Global parameters

% Bayesian observer models

% Noise parameters (~discrimination threshold at given contrast)
params.contrast_sigma = [1 1];
params.contrast_epsilon = [0.05 0.05];

for iParam = 1:numel(params.contrast_sigma)
    params.names{iParam} = 'contrast_sigma';
end
for iParam = 1:numel(params.contrast_epsilon)
    params.names{end+1} = 'contrast_epsilon';
end

params.marginalize_contrasts = true;

% Lapse rate
params.lapse_rate = 0;

% Lapse bias (probability of responding "Left" on a lapse)
params.lapse_bias = 0.5;

% Attentional shift (multiplicative/divisive factor to SIGMA for left/right)
params.attention_factor = 1;

% Softmax parameters for probabilistic mapping from posterior to choice
params.softmax_eta = Inf;     % Deterministic choice, take argmax
params.softmax_bias = 0;      % No bias by default

params.names{end+1} = 'lapse_rate';
extra_features{end+1} = 'lapse';        

%% EXPONENTIAL AVERAGING MODEL

% Function handle to change-point prior
params.runlength_tau = 60;
params.names{end+1} = 'runlength_tau';

if contains(model_name,'nobias')
    params.beta_hyp = [0,0];
else
    % Beta hyperprior on observations
    params.beta_hyp = sqrt([0.1,0.1]);      % Little bias by default

    params.names{end+1} = 'beta_hyp';
    params.names{end+1} = 'beta_hyp';
end

% Assign extra features to model description
if ~isempty(extra_features)
    extra_desc = ' (';
    for iFeat = 1:numel(extra_features)
        extra_desc = [extra_desc, extra_features{iFeat}];
        if iFeat < numel(extra_features)
            extra_desc = [extra_desc,', '];
        end
    end
    extra_desc = [extra_desc,')'];        
    params.model_desc = [params.model_desc,extra_desc];        
end

end