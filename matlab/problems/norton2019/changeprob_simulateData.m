function dataSim = changeprob_simulateData(subID, models, task, simNum)
%CHANGEPROB_SIMULATEDATA Simulates an observer with the same experimental
%parameters of the specified subject(s) using the specified model(s) in the
%specified task(s) numSims times.

% (Detailed documentation to be written.)
%
% Author:   Elyse Norton
% Email:    elyse.norton@gmail.com
% Date:     August 2017

if nargin < 1; subID = []; models = []; task = []; simNum = []; end

models_def = {'fixed', 'idealBayesian', 'exponential', 'RL_probability', ...
    'exponential_conservative', 'RL_probability_conservative', 'RL_criterion', ...
    'subBayesian_rlprior', 'subBayesian_conservative', 'subBayesian_pVec', ...
    'subBayesian_betahyp', 'subBayesian_3param', 'gold', 'gold_nu'}; % Default models

subID_def = {'CWG', 'EGC', 'EHN', 'ERK', 'GK', 'HHL', 'JKT', 'JYZ', 'RND', 'SML', 'SQC'}; % Default subjects in the covert and overt tasks
subID_mixed_def = {'CWG', 'EGC', 'EHN', 'ERK', 'HHL', 'RND', 'SML'}; % Default subjects in the mixed design task

if nargin < 2 || isempty(models)
    models = models_def; % Default simulate all models
end

if nargin < 3 || isempty(task); task = [1 2 3]; end % Default simulate all tasks

if sum(task == 3) == 1 && ~isempty(subID)
    subID_mixed = subID;
end

if isempty(subID)
    subID = subID_def; % Default all subjects for tasks 1 and 2
    subID_mixed = subID_mixed_def; % Default all subjects for task 3
end

if nargin < 4 || isempty(simNum); simNum = 1; end % Default number of simulations is 1/person/task & model

numSims = numel(simNum);

paramBounds_def = [1,30; 1,30; 0,0.1; -Inf,Inf; 0,1; 0,1; 2,200; 0,.5; 0,10; 1.01,5; 1.01,14; 0,1; 0,5];

%% For each subject, model, and task simulate an observer

currentFolder = pwd;

for i = 1:numel(task)
    if task(i) == 1
        currentTask = {'Overt'};
    elseif task(i) == 2
        currentTask = {'Covert'};
    else
        currentTask = {'Mixed'};
    end
    if task(i) ~= 3
        for j = 1:numel(models)
            gridSize = 100; % Default grid size
            idx_model = find(strcmp(models_def, models{j}) == 1);
            if idx_model < 3
                if i == 1
                    parameters = [0 1 0 0 0 0 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 0 0 0 0 0 0 0];
                end
            elseif idx_model == 3 || idx_model == 4 || idx_model == 7
                if i == 1
                    parameters = [0 1 0 0 1 0 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 1 0 0 0 0 0 0 0 0];
                end
            elseif idx_model == 8
                if i == 1
                    parameters = [0 1 0 0 0 0 1 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 1 0 0 0 0 0 0];
                end
            elseif idx_model == 9
                if i == 1
                    parameters = [0 1 0 0 0 1 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 1 0 0 0 0 0 0 0];
                end
            elseif idx_model == 10
                if i == 1
                    parameters = [0 1 0 0 0 0 0 1 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 0 1 0 0 0 0 0];
                end
            elseif idx_model == 11
                if i == 1
                    parameters = [0 1 0 0 0 0 0 0 1 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 0 0 1 0 0 0 0];
                end
            elseif idx_model == 12
                if i == 1
                    parameters = [0 1 0 0 0 0 1 0 1 0 0 0 0];
                else
                    parameters = [1 0 0 0 0 0 1 0 1 0 0 0 0];
                end
                gridSize = 50;
            elseif idx_model == 13
                if i == 1
                    parameters = [0 1 0 0 0 0 0 0 0 1 1 0 0];
                else
                    parameters = [1 0 0 0 0 0 0 0 0 1 1 0 0];
                end
            elseif idx_model == 14
                if i == 1
                    parameters = [0 1 0 0 0 0 0 0 0 1 1 0 1];
                else
                    parameters = [1 0 0 0 0 0 0 0 0 1 1 0 1];
                end
                gridSize = 50;
            else
                if i == 1
                    parameters = [0 1 0 0 1 1 0 0 0 0 0 0 0];
                else
                    parameters = [1 0 0 0 1 1 0 0 0 0 0 0 0];
                end
            end
            NumParams = sum(parameters);
            I_params = find(parameters ~= 0);
            paramBounds = paramBounds_def(I_params,:);
            gridSize = gridSize*ones(1, NumParams);
            for iParam = 1:NumParams
                switch I_params(iParam)
                    case {1, 2}
                        params2fit(iParam,:) = log(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % sigma_noise
                    case 3
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % lambda
                    case 4
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % gamma
                    case 5
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % alpha
                    case 6
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % w
                    case 7
                        params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % Tmax (discrete)
                    case 8
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % Range for minimum probability (pVec(1))
                    case 9
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % beta (hyperprior)
                    case {10, 11}
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % distance between nodes (sqrt of the true deltas)
                    case 12
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % hazard rate
                    case 13
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % nu_p
                end
            end
            for ii = 1:numel(subID)
                idx_subject = find(strcmp(subID_def, subID{ii}) == 1);
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilities_', subID{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb_1:24:17/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                % Choose parameter vector with probability proportional to the joint posterior over parameters for given subject
                postVector = modelPost(:)/sum(modelPost(:));
                % Sort posterior vector
                [postSorted, I_post] = sortrows(postVector);
                % Compute the cumulative sum of sorted vector
                cumPost = cumsum(postSorted);
                for jj = 1:numSims
                    randomSeed = task(i)*idx_model*idx_subject*simNum(jj);
                    rng(randomSeed);
                    % Sample from posterior
                    I_cumsum = find(rand(1,1) < cumPost);
                    I_cumsum = I_cumsum(1); % Choose first element
                    [idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(modelPost), I_post(I_cumsum));
                    for iParam = 1:NumParams
                        simParams(iParam) = params2fit(iParam, idx(iParam)); 
                    end
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || ...
                            strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec') || ...
                            strcmp(models{j}, 'subBayesian_betahyp') || strcmp(models{j}, 'subBayesian_3param')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'gold') || strcmp(models{j}, 'gold_nu')
                        dataSim = changeprob_gold_simulate(data, task(i), simParams);
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, simParams);
                    end
                    dataSim.SimParameters = simParams;
                    dataSim.randomSeed = randomSeed;
                    dataSim.SubjectID = char(subID_def{idx_subject});
                    cd('/Users/elysenorton/Desktop/ChangeProb_1:24:17/ModelSimulations/Simulations_PostSampling');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', ...
                        subID_def{idx_subject}, '_', num2str(simNum(jj)))), 'dataSim');
                    cd(currentFolder);
                    clear simParams
                end
            end
        end
    else
        for j = 1:numel(models)
            gridSize = 100;
            idx_model = find(strcmp(models_def, models{j}) == 1);
            if idx_model < 3
                parameters = [0 1 0 0 0 0 0 0 0 0 0 0 0];
            elseif idx_model == 3 || idx_model == 4 || idx_model == 7
                parameters = [0 1 0 0 1 0 0 0 0 0 0 0 0];
            elseif idx_model == 8
                parameters = [0 1 0 0 0 0 1 0 0 0 0 0 0];
            elseif idx_model == 9
                parameters = [0 1 0 0 0 1 0 0 0 0 0 0 0];
            elseif idx_model == 10
                parameters = [0 1 0 0 0 0 0 1 0 0 0 0 0];
            elseif idx_model == 11
                parameters = [0 1 0 0 0 0 0 0 1 0 0 0 0];
            elseif idx_model == 12
                parameters = [0 1 0 0 0 0 1 0 1 0 0 0 0];
                gridSize = 50;
            elseif idx_model == 13
                parameters = [0 1 0 0 0 0 0 0 0 1 1 0 0];
            elseif idx_model == 14
                parameters = [0 1 0 0 0 0 0 0 0 1 1 0 1];
                gridSize = 50;
            else
                parameters = [0 1 0 0 1 1 0 0 0 0 0 0 0];
            end
            NumParams = sum(parameters);
            I_params = find(parameters ~= 0);
            paramBounds = paramBounds_def(I_params,:);
            gridSize = gridSize*ones(1, NumParams);
            for iParam = 1:NumParams
                switch I_params(iParam)
                    case {1, 2}
                        params2fit(iParam,:) = log(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % sigma_noise
                    case 3
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % lambda
                    case 4
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % gamma
                    case 5
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % alpha
                    case 6
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % w
                    case 7
                        params2fit(iParam,:) = round(linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam))); % Tmax (discrete)
                    case 8
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % Range for minimum probability (pVec(1))
                    case 9
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % beta (hyperprior)
                    case {10, 11}
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % distance between nodes (sqrt of the true deltas)
                    case 12
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % hazard rate
                    case 13
                        params2fit(iParam,:) = linspace(paramBounds(iParam,1), paramBounds(iParam,2), gridSize(iParam)); % nu_p
                end
            end
            for ii = 1:numel(subID_mixed)
                idx_subject = find(strcmp(subID_mixed_def, subID_mixed{ii}) == 1);
                % Load data
                cd('/Users/elysenorton/Desktop/ChangeProb/data');
                load(strcat('ChangingProbabilitiesMixed_', subID_mixed{ii}));
                % Load fit data
                cd('/Users/elysenorton/Desktop/ChangeProb_1:24:17/matlab/ModelFit_data/SMALLlapse');
                load(char(strcat(subID_mixed{ii}, '_', models{j}, '_', currentTask)));
                cd(currentFolder);
                % Choose parameter vector with probability proportional to the joint posterior over parameters for given subject
                postVector = modelPost(:)/sum(modelPost(:));
                % Sort posterior vector
                [postSorted, I_post] = sortrows(postVector);
                % Compute the cumulative sum of sorted vector
                cumPost = cumsum(postSorted);
                for jj = 1:numSims
                    randomSeed = task(i)*idx_model*idx_subject*simNum(jj);
                    rng(randomSeed);
                    % Sample from posterior
                    I_cumsum = find(rand(1,1) < cumPost);
                    I_cumsum = I_cumsum(1); % Choose first element
                    [idx(1),idx(2),idx(3),idx(4),idx(5)] = ind2sub(size(modelPost), I_post(I_cumsum));
                    for iParam = 1:NumParams
                        simParams(iParam) = params2fit(iParam, idx(iParam)); 
                    end
                    if strcmp(models{j}, 'fixed')
                        dataSim = changeprob_fixed_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'idealBayesian') || strcmp(models{j}, 'subBayesian_rlprior') || ...
                            strcmp(models{j}, 'subBayesian_conservative') || strcmp(models{j}, 'subBayesian_pVec') || ...
                            strcmp(models{j}, 'subBayesian_betahyp') || strcmp(models{j}, 'subBayesian_3param')
                        dataSim = changeprob_bocpd_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'exponential') || strcmp(models{j}, 'exponential_conservative')
                        dataSim = changeprob_exp_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'RL_probability') || strcmp(models{j}, 'RL_probability_conservative')
                        dataSim = changeprob_RLprob_simulate(data, task(i), models{j}, simParams);
                    elseif strcmp(models{j}, 'gold') || strcmp(models{j}, 'gold_nu')
                        dataSim = changeprob_gold_simulate(data, task(i), simParams);
                    else
                        dataSim = changeprob_RLcriterion_simulate(data, task(i), models{j}, simParams);
                    end
                    dataSim.SimParameters = simParams;
                    dataSim.randomSeed = randomSeed;
                    dataSim.SubjectID = char(subID_mixed_def{idx_subject});
                    cd('/Users/elysenorton/Desktop/ChangeProb_1:24:17/ModelSimulations/Simulations_PostSampling');
                    save(char(strcat('ChangeProb_Sim_', models{j}, '_', currentTask, '_', subID_mixed_def{idx_subject}, '_', num2str(simNum(jj)))), 'dataSim');
                    cd(currentFolder);
                    clear simParams
                end
            end
        end
    end
end

end

