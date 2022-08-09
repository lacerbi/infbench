function probstruct = infprob_vbmc22(prob,subprob,noise,id,options)

% Problem names (ordered)
problist = {'beta','conbananas','acerbi2012','acerbidokka2018','krajbich2010', ...
    'goris2015b','akrami2018b','wood2010','acerbidokka2018het'};

% Initialize problem structure
if ischar(subprob); D = extractnum(subprob); else; D = subprob; end

switch prob
    case 'beta'
        subprobString = ['P' num2str(D)];
    case {'acerbi2012','acerbidokka2018','krajbich2010','goris2015b','akrami2018b','acerbidokka2018het'}
        subprobString = ['S' num2str(D)];
    otherwise
        subprobString = ['D' num2str(D)];
end

probstruct = initprob(prob,problist,'vbmc22',subprobString,id,D);
% [~,probstruct.Title] = prob;
probstruct.Title = probstruct.Prob;
probstruct.func = ['@(x_,probstruct_) infbench_' probstruct.Prob '(x_(:)'',probstruct_.ProbInfo)'];

probstruct.Noise = noise;
probstruct.NoiseEstimate = 0;       % Function is intrinsically not-noisy

probstruct.ScaleVariables = false;  % Do not rescale
probstruct.PriorType = 'uniform';   % Use uniform prior

switch prob
    case 'beta'
    case 'conbananas'
    case 'wood2010'
    case 'acerbi2012'        
    case 'acerbidokka2018'
    case 'krajbich2010'
    case 'akrami2018b'
    case 'goris2015b'
        probstruct.OptimizeFirst = (D < 100);   % Run optimization        
    case 'acerbidokka2018het'
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
