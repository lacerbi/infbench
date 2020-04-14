function probstruct = infprob_vbmc20(prob,subprob,noise,id,options)

% Problem names (ordered)
problist = {'wood2010','acerbidokka2018','price2018',...
    'krajbich2010','akrami2018','goris2015b','akrami2018b',...
    'nnetcancer'};

% Initialize problem structure
if ischar(subprob); D = extractnum(subprob); else; D = subprob; end

switch prob
    case 'acerbidokka2018'
        subprobString = ['S' num2str(D)];
    case 'wood2010'
        subprobString = ['D' num2str(D)];
    case 'price2018'
        subprobString = ['D' num2str(D)];
    case 'krajbich2010'
        subprobString = ['S' num2str(D)];
    case {'akrami2018','akrami2018b'}
        subprobString = ['S' num2str(D)];
    case 'goris2015b'
        subprobString = ['S' num2str(D)];        
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

probstruct.ScaleVariables = false;  % Do not rescale
probstruct.PriorType = 'uniform';   % Use uniform prior

switch prob
    case 'wood2010'
    case 'acerbidokka2018'
    case {'akrami2018','akrami2018b'}
    case 'goris2015b'
        probstruct.OptimizeFirst = true;    % Run optimization
    case 'nnetcancer'
        probstruct.IntrinsicNoisy = true;
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
