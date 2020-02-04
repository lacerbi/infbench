function probstruct = infprob_vbmc20(prob,subprob,noise,id,options)

% Problem names (ordered)
problist = {'wood2010','nnetcancer'};

% Initialize problem structure
if ischar(subprob); D = extractnum(subprob); else; D = subprob; end

switch prob
    case 'wood2010'
        subprobString = ['D' num2str(D)];
    case 'nnetcancer'
        subprobString = ['D' num2str(D)];        
    otherwise
        subprobString = ['D' num2str(D)];
end

probstruct = initprob(prob,problist,'vbmc20',subprobString,id,D);
% [~,probstruct.Title] = prob;
probstruct.Title = probstruct.Prob;
probstruct.func = ['@(x_,probstruct_) infbench_' probstruct.Prob '(x_(:)'',probstruct_.ProbInfo)'];

probstruct.Noise = noise;
probstruct.NoiseEstimate = 0;       % Function is intrinsically not-noisy

switch prob
    case 'wood2010'
        probstruct.ScaleVariables = false;  % Do not rescale
        probstruct.PriorType = 'uniform';   % Use uniform prior
    case 'nnetcancer'
        probstruct.IntrinsicNoisy = true;
        probstruct.ScaleVariables = false;  % Do not rescale
        probstruct.PriorType = 'uniform';   % Use uniform prior
        probstruct.InferNoise = true;       % Noise is unknown
        probstruct.NoiseEstimate = [];
        probstruct.ComputeTestStatistics = true;
    otherwise
end

%--------------------------------------------------------------------------
function probstruct = initprob(prob,problist,ProbSet,SubProb,id,D)
%INITPROB Initialize problem structure.

if isnumeric(prob); prob = problist{prob}; end
fun = str2func(['infbench_' prob]);
probstruct.ProbSet = ProbSet;
probstruct.Number = find(strcmp(problist,prob),1);
probstruct.Prob = prob;
probstruct.SubProb = SubProb;
probstruct.Id = id;
probstruct.ProbInfo = fun([],[D,id]);
