function [history,post,algoptions] = infalgo_vbmc(algo,algoset,probstruct,getoptions)

if nargin < 4 || isempty(getoptions)
    getoptions = false;
end

algoptions = vbmc('all');                   % Get default settings

ControlRunFlag = false;     % Do NOT run in control mode

if nargin < 3 || isempty(probstruct)
    probstruct.MinFunEvals = algoptions.MinFunEvals;
    probstruct.MaxFunEvals = algoptions.MaxFunEvals;
    probstruct.D = 5;
    probstruct.Debug = false;
    probstruct.InitPoint = zeros(1,probstruct.D);
    probstruct.Noise = [];
    probstruct.IntrinsicNoisy = false;
else
    algoptions.MinFunEvals = probstruct.MaxFunEvals;
    algoptions.MaxFunEvals = probstruct.MaxFunEvals;
end

% Default VBMC options as of 1.0.9
algoptions.FunEvalsPerIter = 5;
algoptions.AcqFcn = '@vbmc_acqskl';
algoptions.SearchAcqFcn = '@acqf_vbmc';
algoptions.SearchOptimizer = 'cmaes';
algoptions.MinIter = 0;     % No limits on iterations
algoptions.MaxIter = Inf;
algoptions.MinFinalComponents = 0;
algoptions.WarpNonlinear = 'off';   % No nonlinear warping for now
algoptions.BestSafeSD = 5;
algoptions.BestFracBack = 0.25;
algoptions.Diagnostics = 'on';
algoptions.InitDesign = 'plausible';    % Initial design uniform in plausible box
algoptions.EmpiricalGPPrior = 'yes';
algoptions.TolStableWarmup = 15;
algoptions.TolStableExceptions = 1/8;
algoptions.RankCriterion = true;
algoptions.EmpiricalGPPrior = true;
algoptions.gpQuadraticMeanBound = false;
algoptions.WarmupOptions = [];
algoptions.NSentFastBoost = [];
algoptions.NSentFineBoost = [];
algoptions.ActiveVariationalSamples = 0;
algoptions.GPTrainNinit = 1024;
algoptions.DetEntropyAlpha = 0;
algoptions.VariationalInitRepo = false;
algoptions.UpperGPLengthFactor = 0;
algoptions.TolGPVarMCMC = 1e-4;
algoptions.PosteriorMCMC = 0;
algoptions.VarThresh = Inf;
algoptions.ActiveImportanceSamplingMCMCSamples = 100;
algoptions.MaxRepeatedObservations = 0;
algoptions.WarpMinK = 5;
algoptions.ActiveSampleVPUpdate = false;
algoptions.ActiveSampleGPUpdate = false;
algoptions.ActiveSampleFullUpdatePastWarmup = 2;
algoptions.ActiveSampleFullUpdateThreshold = 3;
algoptions.ActiveSamplefESSThresh = 1;
algoptions.NoiseShaping = false;
algoptions.NoiseShapingThreshold = 10*probstruct.D;
algoptions.NoiseShapingFactor = 0.05;
algoptions.WarpTolImprovement = 0.1;
algoptions.WarpUndoCheck = true;
algoptions.PruningThresholdMultiplier = @(K) 1/sqrt(K);
algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint);
algoptions.TolStableExcptFrac = 0.2;
algoptions.TolStableCount = 60;
algoptions.WarmupCheckMax = true;
algoptions.SGDStepSize = 0.005;
algoptions.NSentFine = '@(K) 2^12*K';
algoptions.NSentFast = 0;
algoptions.gpMeanFun = 'negquad';
algoptions.GPTrainInitMethod = 'rand';
algoptions.GPTrainNinitFinal = 64;
algoptions.StopWarmupReliability = 100;
algoptions.WarmupKeepThresholdFalseAlarm = '100*(D+2)';
algoptions.NSentActive = '@(K) 20*K.^(2/3)';
algoptions.NSentBoost = '@(K) 200*K.^(2/3)';
algoptions.SkipActiveSamplingAfterWarmup = 0;
algoptions.TolGPNoise = sqrt(1e-5);
algoptions.WarpRotoCorrThresh = 0.05;
algoptions.EmpiricalGPPrior = 0;
algoptions.GPLengthPriorMean = 'sqrt(D/6)';
algoptions.GPLengthPriorStd = 0.5*log(1e3);
algoptions.ActiveSearchBound = 2; 
algoptions.BoxSearchFrac = 0.25;
algoptions.GPTolOpt = 1e-5;
algoptions.GPTolOptMCMC = 1e-2;
algoptions.gpQuadraticMeanBound = 1;
algoptions.NSent = '@(K) 100*K.^(2/3)';
algoptions.MaxIterStochastic = '100*(2+nvars)';
algoptions.SearchMaxFunEvals = '500*(D+2)';
algoptions.GPSampleThin = 5;
algoptions.StableGPvpK = Inf;
algoptions.WarmupKeepThreshold = '10*D';
algoptions.RecomputeLCBmax = true;
algoptions.MinFinalComponents = 50;
algoptions.WarpRotoScaling = true;
algoptions.BoundedTransform = 'logit';

if probstruct.Debug
    algoptions.TrueMean = probstruct.Post.Mean;
    algoptions.TrueCov = probstruct.Post.Cov;
end

% Use prior as proposal function
algoptions.ProposalFcn = @(X_) exp(infbench_lnprior(X_,probstruct));

% Default options for variational active sampling
if numel(algoset) >= 3 && strcmpi(algoset(1:3),'vas')
    algoptions.VarActiveSample = true;
    algoptions.VariableMeans = false;
    algoptions.Plot = 'on';
    algoptions.NSent = 0;
    algoptions.NSentFast = 0;
    algoptions.NSentFine = '@(K) 2^15*round(sqrt(K))';
    algoptions.DetEntTolOpt = 1e-3;
    algoptions.EntropySwitch = true;
    algoptions.DetEntropyMinD = 0;
    algoptions.EntropyForceSwitch = Inf;
    algoptions.TolWeight = 0;
    algoptions.NSelbo = '@(K) 50*sqrt(K)';
    algoptions.SearchAcqFcn = @acqvasreg_vbmc;
    algoptions.NSsearch = 2^7;
    algoptions.Warmup = false;
    algoptions.FunEvalsPerIter = 1;
end

% Changing these options reduces performance on noisy datasets (for
% uncertainty sampling acquisition function)
% algoptions.StableGPvpK = 10;
% algoptions.GPSampleThin = 3;
% algoptions.WarmupKeepThreshold = '100*(D+2)';

% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = true; algoptions.Plot = 'on'; algoptions.FeatureTest = true;
    case {1,'base'}; algoset = 'base';                                                      % Use defaults
    case {2,'acqusreg'}; algoset = 'acqusreg'; algoptions.SearchAcqFcn = @acqus_vbmc;       % Vanilla uncertainty sampling
    case {3,'acqproreg'}; algoset = 'acqproreg'; algoptions.SearchAcqFcn = @acqf_vbmc;      % Prospective uncertainty sampling
    case {4,'control'}; algoset = 'control'; ControlRunFlag = true;                         % Control experiment
    case {5,'test'}; algoset = 'test'; algoptions.FeatureTest = true;                       % Test feature
    case {6,'narrow'}; algoset = 'narrow'; algoptions.InitDesign = 'narrow';                % Narrow initialization
    case {7,'control2'}; algoset = 'control2'; ControlRunFlag = true;                       % Control experiment, repeated
    case {8,'test2'}; algoset = 'test2'; algoptions.FeatureTest = true;                     % Test feature (second case)
    case {9,'fastdebug'}; algoset = 'fastdebug'; algoptions.Debug = true; algoptions.Plot = 'off'; algoptions.FeatureTest = true; algoptions.MinFinalComponents = 0;

    % Fixed number of mixture components
    case {11,'K1'}; algoset = 'K1'; algoptions.Kfun = 1; algoptions.KfunMax = 1; algoptions.Kwarmup = 1;
    case {12,'K2'}; algoset = 'K2'; algoptions.Kfun = 2; algoptions.KfunMax = 2; algoptions.Kwarmup = 2;
    case {13,'K5'}; algoset = 'K5'; algoptions.Kfun = 5; algoptions.KfunMax = 5; algoptions.Kwarmup = 5;
        
    % Ongoing research and testing
    case {21,'acqfregvlnn'}; algoset = 'acqfregvlnn'; algoptions.SearchAcqFcn = @vbmc_acqfregvlnn;
    case {22,'acqfregvsqrtn'}; algoset = 'acqfregvsqrtn'; algoptions.SearchAcqFcn = @vbmc_acqfregvsqrtn;
    case {23,'acqfregt'}; algoset = 'acqfregt'; algoptions.SearchAcqFcn = @vbmc_acqfregt;   % annealed prospective uncertainty sampling
    case {24,'se'}; algoset = 'se'; algoptions.gpMeanFun = 'se';           
    case {25,'const'}; algoset = 'const'; algoptions.gpMeanFun = 'const';        
    case {26,'acqf2reg'}; algoset = 'acqf2reg'; algoptions.SearchAcqFcn = @vbmc_acqf2reg;
    case {27,'acqpropnew'}; algoset = 'acqpropnew'; algoptions.SearchAcqFcn = @acqpropreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {28,'bowarmup'}; algoset = 'bowarmup';
        w.SearchCacheFrac = 0.1; w.HPDSearchFrac = 0.9; w.HeavyTailSearchFrac = 0; w.MVNSearchFrac = 0; w.SearchAcqFcn = @acqpropreg_vbmc; w.StopWarmupThresh = 0.1; w.SearchCMAESVPInit = false;
        algoptions.WarmupOptions = w; algoptions.Plot = 0;
        algoptions.TolStableWarmup = 55; algoptions.BOWarmup = true; algoptions.NSgpMaxWarmup = 8;
        algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 220 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {29,'acqfheavy'}; algoset = 'acqfheavy'; algoptions.Plot = 1; algoptions.SearchAcqFcn = @acqfheavyreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {30,'gpswitch'}; algoset = 'gpswitch'; algoptions.Plot = 0; algoptions.StableGPSampling = '100+10*nvars'; algoptions.WarmupOptions.StableGPSampling = '200+10*nvars'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {31,'noswitch'}; algoset = 'noswitch'; algoptions.Plot = 0; algoptions.StableGPSampling = 'Inf'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {32,'hiK'}; algoset = 'hiK'; algoptions.Plot = 0; algoptions.KfunMax = '@(N) min(N,100)'; algoptions.AdaptiveK = 5; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {33,'outwarp'}; algoset = 'outwarp'; algoptions.Plot = 0; algoptions.gpOutwarpFun = 'outwarp_negscaledpow'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {34,'outwarp2'}; algoset = 'outwarp2'; algoptions.Plot = 0; algoptions.gpOutwarpFun = 'outwarp_negpow'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {35,'outwarp3'}; algoset = 'outwarp3'; algoptions.Plot = 0; algoptions.gpOutwarpFun = 'outwarp_negpow'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {36,'outwarp4'}; algoset = 'outwarp4'; algoptions.Plot = 0; algoptions.gpOutwarpFun = 'outwarp_negpowc1'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {37,'outwarp5'}; algoset = 'outwarp5'; algoptions.Plot = 0; algoptions.gpOutwarpFun = 'outwarp_negpowc1'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {38,'outwarp6'}; algoset = 'outwarp6'; algoptions.Plot = 0; algoptions.WarmupKeepThreshold = 50*numel(probstruct.InitPoint); algoptions.gpOutwarpFun = 'outwarp_negpowc1'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {39,'T2'}; algoset = 'T2'; algoptions.Temperature = 2; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {40,'T2tolw'}; algoset = 'T2tolw'; algoptions.TolWeight = 1e-3; algoptions.Temperature = 2; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {41,'tolw'}; algoset = 'tolw'; algoptions.TolWeight = 1e-3; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {42,'T2tolwwarp'}; algoset = 'T2tolw'; algoptions.TolWeight = 1e-3; algoptions.gpOutwarpFun = 'outwarp_negpowc1';  algoptions.Temperature = 2; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {43,'acqvttolw'}; algoset = 'acqvttolw'; algoptions.Plot = 0; algoptions.SearchAcqFcn = @acqfregvrnd_vbmc; algoptions.TolWeight = 1e-3; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {44,'T2tolwacqv'}; algoset = 'T2tolwacqv'; algoptions.TolWeight = 1e-3; algoptions.SearchAcqFcn = @acqfregvrnd_vbmc; algoptions.Temperature = 2; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {45,'acqiqrtolw'}; algoset = 'acqiqrtolw'; algoptions.Plot = 0; algoptions.TolWeight = 1e-3; algoptions.SearchAcqFcn = @acqiqrreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {46,'gp2'}; algoset = 'gp2'; algoptions.Plot = 1; algoptions.SeparateSearchGP = 1; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {47,'noiseshaping'}; algoset = 'noiseshaping'; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '20*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {48,'acqim'}; algoset = 'acqim'; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {49,'noiseshaping2'}; algoset = 'noiseshaping2'; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {50,'noiseshaping2im'}; algoset = 'noiseshaping2im'; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {51,'acq2'}; algoset = 'acq2'; algoptions.SearchAcqFcn = {@acqeig_vbmc,@acqf_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {54,'acq2shaping'}; algoset = 'acq2shaping'; algoptions.SearchAcqFcn = {@acqf_vbmc,@acqeig_vbmc}; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {55,'acq2hedge'}; algoset = 'acq2hedge'; algoptions.AcqHedge = 1; algoptions.SearchAcqFcn = {@acqf_vbmc,@acqmireg_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {56,'acq2shapinghedge'}; algoset = 'acq2shapinghedge'; algoptions.AcqHedge = 1; algoptions.SearchAcqFcn = {@acqf_vbmc,@acqmireg_vbmc}; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {57,'acq2new'}; algoset = 'acq2new'; algoptions.SearchAcqFcn = {@acqmireg_vbmc,@acqf_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {60,'step1migps'}; algoset = 'acqstep1migps'; algoptions.FunEvalsPerIter = 1; algoptions.FunEvalStart = 'D'; algoptions.KfunMax = @(N) N; algoptions.SeparateSearchGP = 1; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {61,'actfull'}; algoset = 'actfull'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = true;
    case {62,'up4'}; algoset = 'up4'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = true; algoptions.NSent = '@(K) 100*K.^(2/3)'; algoptions.NSentBoost = '@(K) 200*K.^(2/3)'; algoptions.SkipActiveSamplingAfterWarmup = 0; algoptions.StableGPvpK = 10;
    case {63,'fast'}; algoset = 'fast'; algoptions = newdefaults; algoptions.NSent = '@(K) 100*K.^(2/3)'; algoptions.NSentBoost = '@(K) 200*K.^(2/3)'; algoptions.SkipActiveSamplingAfterWarmup = 0; algoptions.StableGPvpK = 10;
    case {64,'actfull2'}; algoset = 'actfull2'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = true; algoptions.GPTolOptActive = 1e-2;
    case {65,'postmcmc'}; algoptions = newdefaults; algoptions.PosteriorMCMC = 1e3;
    case {66,'trust'}; algoset = 'trust'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfregtr_vbmc;
    case {67,'trust2'}; algoset = 'trust2'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmaxiqrregtr_vbmc; algoptions.Plot = 1;
    case {68,'new'}; algoset = 'new'; algoptions = newdefaults; algoptions.TolGPNoise = sqrt(1e-5); algoptions.GPLengthPriorMean = 'sqrt(D/6)'; algoptions.GPLengthPriorStd = 0.5*log(1e3); algoptions.ActiveSearchBound = 2; algoptions.BoxSearchFrac = 0.25;
    case {69,'new2'}; algoset = 'new'; algoptions = newdefaults; algoptions.TolGPNoise = sqrt(1e-5); algoptions.GPLengthPriorMean = 'sqrt(D/6)'; algoptions.GPLengthPriorStd = 0.5*log(1e3); algoptions.ActiveSearchBound = 2; algoptions.BoxSearchFrac = 0.25; algoptions.gpMeanFun = 'negquad';
    case {70,'new3'}; algoset = 'new'; algoptions = newdefaults; algoptions.TolGPNoise = sqrt(1e-4); algoptions.GPLengthPriorMean = 'sqrt(D/6)'; algoptions.GPLengthPriorStd = 0.5*log(1e3); algoptions.ActiveSearchBound = 2; algoptions.BoxSearchFrac = 0.25; algoptions.gpMeanFun = 'negquad'; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveImportanceSamplingMCMCSamples = 200;
    case {71,'new4'}; algoset = 'new'; algoptions = newdefaults; algoptions.TolGPNoise = sqrt(1e-4); algoptions.GPLengthPriorMean = 'sqrt(D/6)'; algoptions.GPLengthPriorStd = 0.5*log(1e3); algoptions.ActiveSearchBound = 2; algoptions.BoxSearchFrac = 0.25; algoptions.gpMeanFun = 'negquad'; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0;
    case {72,'new5'}; algoset = 'new'; algoptions = newdefaults; algoptions.TolGPNoise = sqrt(1e-4); algoptions.GPLengthPriorMean = 'sqrt(D/6)'; algoptions.GPLengthPriorStd = 0.5*log(1e3); algoptions.ActiveSearchBound = 2; algoptions.BoxSearchFrac = 0.25; algoptions.gpMeanFun = 'negquad'; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveImportanceSamplingMCMCSamples = 200; 
        algoptions.WarmupOptions.NSgpMaxWarmup = 8; algoptions.WarmupOptions.SearchAcqFcn = @acqf_vbmc;
    case {73,'warp'}; algoset = 'warp'; algoptions = newdefaults; algoptions.WarpRotoScaling = 1; algoptions.WarpNonlinear = 1;
    case {74,'roto'}; algoset = 'roto'; algoptions = newdefaults; algoptions.WarpRotoScaling = 1;
    case {75,'negquadse'}; algoset = 'negquadse'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadse';
    case {76,'quadmix'}; algoset = 'quadmix'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadmix';
    case {80,'rotoup'}; algoset = 'rotoup'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1;
    case {81,'acqbothup'}; algoset = 'acqbothup'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1; algoptions.WarmupOptions.ActiveSampleFullUpdate = 0; algoptions.WarmupOptions.SearchAcqFcn = '@acqf_vbmc';
    case {82,'acqboth'}; algoset = 'acqboth'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxMain = 0; algoptions.WarpRotoScaling = 1; algoptions.WarmupOptions.SearchAcqFcn = '@acqf_vbmc';
    
    % Final testing
    case {101,'oldsettings'}; algoset = 'oldsettings';
    case {102,'newbase'}; algoset = 'newbase'; algoptions = newdefaults; algoptions.SearchAcqFcn = '@acqf_vbmc'; algoptions.WarpRotoScaling = 0; algoptions.MinFinalComponents = 0;
    case {103,'renewdef'}; algoset = 'renewdef'; algoptions = algoptions;
    case {104,'renewdefnoroto'}; algoset = 'renewdef'; algoptions = algoptions; algoptions.WarpRotoScaling = false;
        
    % New defaults
    case {150,'newdef'}; algoset = 'newdef'; algoptions = newdefaults;
    case {151,'newdef2'}; algoset = 'newdef2'; algoptions = newdefaults;    % Current best
    case {152,'newdef3'}; algoset = 'newdef3'; algoptions = newdefaults;
    case {154,'newdef4'}; algoset = 'newdef4'; algoptions = newdefaults; algoptions.Plot = 1;
    case {160,'newdefdebug'}; algoset = 'newdefdebug'; algoptions = newdefaults; algoptions.MinFinalComponents = 0;
        
    % Noisy paper (base and variants)
    case {201,'eig'}; algoset = 'eig';      algoptions.SearchAcqFcn = @acqeig_vbmc;             algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {202,'npro'}; algoset = 'npro';    algoptions.SearchAcqFcn = @acqfsn2_vbmc;            algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {203,'viqr'}; algoset = 'viqr';    algoptions.SearchAcqFcn = @acqviqr_vbmc;            algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 100;
    case {204,'imiqr'}; algoset = 'imiqr';  algoptions.SearchAcqFcn = @acqimiqr_vbmc;           algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 100; algoptions.ActiveImportanceSamplingMCMCThin = 5;
    case {211,'viqrnoroto'}; algoset = 'viqrnoroto'; algoptions.SearchAcqFcn = @acqviqr_vbmc;   algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 0; algoptions.MinFinalComponents = 0;
    case {212,'viqrnogp'}; algoset = 'viqrnogp'; algoptions.SearchAcqFcn = @acqviqr_vbmc;       algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0;

    % Old names
    case {251,'renewdefmipluswup4gpsvp'}; algoset = 'renewdefmipluswup4gpsvp';          algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc;          algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {260,'renewdefimiqrpluswup5noacq'}; algoset = 'renewdefimiqrpluswup5noacq';    algoptions = algoptions; algoptions.SearchAcqFcn = @acqfsn2_vbmc;        algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {259,'renewdefvarimiqrpluswup5fast'}; algoset = 'renewdefvarimiqrpluswup5fast';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 100;
    case {264,'renewdefimiqrplus5longvpgps'}; algoset = 'renewdefimiqrplus5longvpgps';  algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc;       algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 100; algoptions.ActiveImportanceSamplingMCMCThin = 5;
        
    % Noise
    case {231,'acqsn2'}; algoset = 'acqsn2'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfsn2reg_vbmc; algoptions.Plot = 0; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.WarmupKeepThreshold = '1e3*(nvars+2)'; algoptions.WarmupKeepThresholdFalseAlarm = '1e3*(nvars+2)'; algoptions.MaxRepeatedObservations = 0; algoptions.PosteriorMCMC = 2e4; algoptions.VarThresh = 1;
    % case {202,'heur'}; algoset = 'heur'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfsn2regtrheur_vbmc; algoptions.Plot = 1; algoptions.ActiveSampleFullUpdate = 1;
    case {235,'acqnoise'}; algoset = 'acqnoise'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfsn2_vbmc; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1;
    case {237,'acqnoisemcmc'}; algoset = 'acqnoisemcmc'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfsn2_vbmc; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1; algoptions.PosteriorMCMC = 2e4;    
    case {238,'renewdefnoise'}; algoset = 'renewdefnoise'; algoptions = algoptions; algoptions.SearchAcqFcn = '@acqfsn2_vbmc'; algoptions.WarpRotoScaling = 0; algoptions.MinFinalComponents = 0;
    case {239,'renewdefnoiseplus'}; algoset = 'renewdefnoiseplus'; algoptions = algoptions; algoptions.SearchAcqFcn = '@acqfsn2_vbmc'; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveSampleFullUpdate = 2;

    case {240,'renewdefimiqrplus'}; algoset = 'renewdefimiqrplus'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {241,'renewdefimiqrpluswup'}; algoset = 'renewdefimiqrpluswup'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.WarmupOptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;        
    case {248,'renewdefmipluswup4'}; algoset = 'renewdefmipluswup4'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {249,'renewdefmipluswup4gps'}; algoset = 'renewdefmipluswup4gps'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {250,'renewdefmipluswup4vp'}; algoset = 'renewdefmipluswup4'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {254,'renewdefmipluswup4gpsvplcbmix'}; algoset = 'renewdefmipluswup4gpsvp'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc; algoptions.WarmupOptions.SearchAcqFcn = @acqfsn2_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.RecomputeLCBmax = true;
    case {256,'renewdefimiqrpluswup5'}; algoset = 'renewdefimiqrpluswup5'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {257,'renewdefimiqrvppluswup5'}; algoset = 'renewdefimiqrvppluswup5'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqrvp_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {258,'renewdefvarimiqrpluswup5vp'}; algoset = 'renewdefvarimiqrpluswup5'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
        
    case {261,'renewdefvarimiqrpluswup5vpfmincon'}; algoset = 'renewdefvarimiqrpluswup5'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.NSsearch = 2^11;
    case {262,'renewdefimiqrplus5long'}; algoset = 'renewdefimiqrplus5long'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 500; algoptions.ActiveImportanceSamplingMCMCThin = 11;
    case {263,'renewdefimiqrplus5longvp'}; algoset = 'renewdefimiqrplus5long'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 500; algoptions.ActiveImportanceSamplingMCMCThin = 11;
    case {265,'renewdeffimiqrpluswup5vp'}; algoset = 'renewdeffimiqrpluswup5vp'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqfimiqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCThin = 5;
    case {266,'renewdeffimiqrpluswup5vpfess'}; algoset = 'renewdeffimiqrpluswup5vpfess'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqfimiqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCThin = 5; algoptions.ActiveSamplefESSThresh = 0.75;
    case {267,'renewdefvarimiqrmixpluswup5vp'}; algoset = 'renewdefvarimiqrmixpluswup5vp'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqvarimiqrmix_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
    case {268,'renewdefvarimiqrpluswup5fastfess'}; algoset = 'renewdefvarimiqrpluswup5fastfess'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveSamplefESSThresh = 0.75;

    case {270,'viqrnorminv'}; algoset = 'viqrnorminv';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.BoundedTransform = 'norminv';
    case {271,'viqrstudent4'}; algoset = 'viqrstudent4';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.BoundedTransform = 'student4';
    case {272,'viqr500'}; algoset = 'viqr500';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 500;
    case {274,'viqrlaplace'}; algoset = 'viqrlaplace'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.GPHypSampler = 'laplace';
    case {275,'viqrnpv'}; algoset = 'viqrnpv'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.GPHypSampler = 'npv';
    case {277,'viqrnogplate'}; algoset = 'viqrnogplate'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.NSgpMaxMain = 0;
    case {278,'viqr2gp'}; algoset = 'viqr2gp'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.DoubleGP = true;
    case {279,'viqrnegquadmix'}; algoset = 'renewdefvarimiqrpluswup5fast';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.gpMeanFun = 'negquadmix';
    case {280,'viqrnotrim'}; algoset = 'viqrnotrim';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.WarmupKeepThresholdFalseAlarm = Inf; algoptions.WarmupKeepThreshold = Inf;        
    case {281,'eignotrim'}; algoset = 'eignotrim';algoptions = algoptions; algoptions.SearchAcqFcn = @acqeig_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.WarmupKeepThresholdFalseAlarm = Inf; algoptions.WarmupKeepThreshold = Inf;        
    case {282,'npronotrim'}; algoset = 'npronotrim';algoptions = algoptions; algoptions.SearchAcqFcn = @acqfsn2_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.WarmupKeepThresholdFalseAlarm = Inf; algoptions.WarmupKeepThreshold = Inf;        
    case {283,'imiqrnotrim'}; algoset = 'imiqrnotrim';algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.WarmupKeepThresholdFalseAlarm = Inf; algoptions.WarmupKeepThreshold = Inf; algoptions.ActiveImportanceSamplingMCMCThin = 5;
    case {284,'viqrnotrimshape'}; algoset = 'viqrnotrimshape';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.WarmupKeepThresholdFalseAlarm = Inf; algoptions.WarmupKeepThreshold = Inf; algoptions.NoiseShaping = 1;
    case {285,'viqrt2'}; algoset = 'viqrt2';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.Temperature = 2;
    case {286,'viqrelcbo'}; algoset = 'viqrelcbo';algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc;    algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.OptimisticVariationalBound = 0.6745; algoptions.ELCBOWeight = -0.6745;

    % Information-theoretic
    case {302,'acqmi'}; algoset = 'acqmi'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; % Needs to be rerun on base/no-noise
    case {303,'acqmidebug'}; algoset = 'acqmidebug'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveVariationalSamples = 100; algoptions.ActiveSampleFullUpdate = 1; algoptions.Plot = 1;
    case {310,'acqsn'}; algoset = 'acqsn'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqfsn2_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.SampleExtraVPMeans = '@(K)10+K'; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.PosteriorMCMC = 2e4; algoptions.Plot = 1; algoptions.WarmupKeepThreshold = '50*(nvars+2)'; algoptions.WarmupKeepThresholdFalseAlarm = '100*(nvars+2)';
%    case {311,'acqopt1'}; algoset = 'acqopt1'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.SampleExtraVPMeans = '@(K)10+K'; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.PosteriorMCMC = 2e4; algoptions.Plot = 1; algoptions.WarmupKeepThreshold = '50*(nvars+2)'; algoptions.WarmupKeepThresholdFalseAlarm = '100*(nvars+2)'; algoptions.OptimisticVariationalBound = 1;
    case {331,'acqimi2'}; algoset = 'acqmi'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.ActiveSampleFullUpdate = 0; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.Plot = 0; algoptions.WarmupKeepThreshold = '1e3*(nvars+2)'; algoptions.WarmupKeepThresholdFalseAlarm = '1e3*(nvars+2)'; algoptions.MaxRepeatedObservations = 0; algoptions.PosteriorMCMC = 2e4;
    case {332,'acqmaxiqr'}; algoset = 'acqmi'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmaxiqrreg_vbmc; algoptions.ActiveSampleFullUpdate = 2; algoptions.Plot = 1; algoptions.WarmupKeepThreshold = 'Inf'; algoptions.WarmupKeepThresholdFalseAlarm = 'Inf'; algoptions.MaxRepeatedObservations = 0; algoptions.PosteriorMCMC = 2e4;
    case {333,'acqimi3'}; algoset = 'acqimi3'; algoptions = newdefaults; algoptions.FunEvalStart = 10; algoptions.Plot = 0; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.Plot = 0; algoptions.WarmupKeepThreshold = '1e3*(nvars+2)'; algoptions.WarmupKeepThresholdFalseAlarm = '1e3*(nvars+2)'; algoptions.MaxRepeatedObservations = 0; algoptions.PosteriorMCMC = 2e4; algoptions.VarThresh = 1;
           
%    case {350,'acqmivar'}; algoset = 'acqmivar'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.Plot = 1; ...
%            algoptions.VariableMeans = 0; algoptions.NSent = 0; algoptions.NSentActive = 0; algoptions.NSentBoost = 0; algoptions.NSentFine = '@(K) 200*K.^(2/3)'; algoptions.NSentFineActive = '@(K) 200*K.^(2/3)'; algoptions.Warmup = 0;
    
    % Integrated mean function (tried a lot of variants - all unstable)
%    case {400,'intmean'}; algoset = 'intmean'; algoptions = newdefaults; algoptions.gpIntMeanFun = 3; algoptions.gpMeanFun = 'zero'; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.Plot = 0;

    % Extra stuff
    case {401,'t2'}; algoset = 't2'; algoptions = algoptions; algoptions.Temperature = 2;
    case {402,'t2viqr'}; algoset = 't2viqr'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.ActiveImportanceSamplingMCMCSamples = 100; algoptions.Temperature = 2;

    % Tests
    case {501,'imiqrplusfit'}; algoset = 'imiqrplusfit'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqimiqr_vbmc; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.ActiveSampleFullUpdate = 2; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50; algoptions.FitnessShaping = true; algoptions.OutwarpThreshBase = '20*(nvars+1)';

    % Tests 2022
    case {600,'basenew'}; algoset = 'basenew';
    case {601,'probit'}; algoset = 'probit'; algoptions.BoundedTransform = 'probit';
    case {602,'noisynew'}; algoset = 'noisynew'; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true;
    case {603,'noisyprobit'}; algoset = 'noisyprobit'; algoptions.BoundedTransform = 'probit'; algoptions.SearchAcqFcn = @acqviqr_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true;
             
    % Variational active sampling
    case {1000,'vas'}; algoset = 'vas'; 
    
    % New features    
    case {10000,'newfeat'}; algoset = 'newfeat'; algoptions = algoptions; algoptions.SearchAcqFcn = @acqbald_vbmc; algoptions.ActiveSampleGPUpdate = true; algoptions.ActiveSampleVPUpdate = true; algoptions.WarpRotoScaling = 1; algoptions.MinFinalComponents = 50;
        
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

% Increase base noise with noisy functions
if ~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy
    algoptions.UncertaintyHandling = 'on';
    
    if probstruct.InferNoise
        algoptions.SpecifyTargetNoise = false;        
        NoiseEstimate = probstruct.NoiseEstimate;
        % if isempty(NoiseEstimate); NoiseEstimate = 1; end
        if ~isempty(NoiseEstimate)
            algoptions.NoiseSize = NoiseEstimate(1);
        else
            algoptions.NoiseSize = [];            
        end
    else
        algoptions.SpecifyTargetNoise = true;
    end    
    
    algoptions.TolStableCount = ceil(algoptions.TolStableCount*1.5);
    % algoptions.TolStableWarmup = algoptions.TolStableWarmup*2;    
else
    algoptions.UncertaintyHandling = 'off';
end

% If GETOPTIONS is true, just return the algorithm options
if getoptions
    history = algoptions;
    post = probstruct;
    return;
end

PLB = probstruct.PLB;
PUB = probstruct.PUB;
LB = probstruct.LB;
UB = probstruct.UB;
x0 = probstruct.InitPoint;
D = size(x0,2);

% Add log prior to function evaluation 
% (the current version of VBMC is agnostic of the prior)
probstruct.AddLogPrior = true;

algo_timer = tic;
MinFinalComponents = algoptions.MinFinalComponents; % Skip final boosting here, do it later
algoptions.MinFinalComponents = 0;
[vp,elbo,elbo_sd,exitflag,~,~,output,stats] = ...
    vbmc(@(x) infbench_func(x,probstruct),x0,LB,UB,PLB,PUB,algoptions);
TotalTime = toc(algo_timer);
gp = vp.gp;
algoptions.MinFinalComponents = MinFinalComponents;

% Get preprocessed OPTIONS struct
options_vbmc = setupoptions_vbmc(D,algoptions,algoptions);

if ~ControlRunFlag
    
    history = infbench_func(); % Retrieve history
    history.scratch.output = output;
    history.TotalTime = TotalTime;
    
    % Store computation results (ignore points discarded after warmup)
    history.Output.X = output.X_orig(output.X_flag,:);
    history.Output.y = output.y_orig(output.X_flag);
    post.lnZ = elbo;
    post.lnZ_var = elbo_sd^2;
    fprintf('Calculating VBMC output at iteration...\n');
    fprintf('%d..',0);
    [post.gsKL,post.Mean,post.Cov,post.Mode,post.MTV,post.samples,post.Test] = ...
        computeStats(vp,gp,probstruct,algoptions.PosteriorMCMC,algoptions.VarThresh);

    % Return estimate, SD of the estimate, and gauss-sKL with true moments
    Nticks = numel(history.SaveTicks);
    for iIter = 1:Nticks
        fprintf('%d..',iIter);
        
        idx = find(stats.funccount == history.SaveTicks(iIter),1);
        if isempty(idx); continue; end
        
        history.ElapsedTime(iIter) = stats.timer(idx).totalruntime;        
        
        % Compute variational solution that VBMC would return at a given iteration
        t = tic();
        [vp,~,~,idx_best] = ...
            best_vbmc(stats,idx,algoptions.BestSafeSD,algoptions.BestFracBack,algoptions.RankCriterion,0);
        [vp,elbo,elbo_sd] = ...
            finalboost_vbmc(vp,idx_best,[],stats,options_vbmc);
        gp = stats.gp(idx_best);    % Get GP corresponding to chosen iter
        
        % Convert training variational posterior to real variational posterior
        vp = vptrain2real(vp,1);
        elbo = vp.stats.elbo;
        elbo_sd = vp.stats.elbo_sd;
        
        % Take actual runtime at the end of each iteration and add boost time
        history.ElapsedTime(iIter) = stats.timer(idx).totalruntime + toc(t);
        
        history.Output.N(iIter) = history.SaveTicks(iIter);
        history.Output.lnZs(iIter) = elbo;
        history.Output.lnZs_var(iIter) = elbo_sd^2;
        [gsKL,Mean,Cov,Mode,MTV,~,Test] = computeStats(vp,gp,probstruct,algoptions.PosteriorMCMC,algoptions.VarThresh);
        history.Output.Mean(iIter,:) = Mean;
        history.Output.Cov(iIter,:,:) = Cov;
        history.Output.gsKL(iIter) = gsKL;
        history.Output.Mode(iIter,:) = Mode;
        history.Output.MTV(iIter,:) = MTV;
        history.Output.Test{iIter} = Test;
    end
    fprintf('\n');
else
    % Control condition -- run VBMC as normal but compute marginal likelihood
    % and posterior via other methods

    % Add WSABI algorithm to MATLAB path
    BaseFolder = fileparts(mfilename('fullpath'));
    AlgoFolder = 'wsabi';
    addpath(genpath([BaseFolder filesep() AlgoFolder]));
        
    history = infbench_func(); % Retrieve history
    
    % Start warmup
    
    % Store all points (no warmup pruning)
    X = output.X_orig(1:output.Xn,:);
    y = output.y_orig(1:output.Xn);
    Niter = find(size(X,1) == history.SaveTicks,1);
    N = history.SaveTicks(1:Niter);
    
    % Find when warm-up ends
    idx = find(stats.warmup == 0,1);
    if isempty(idx); endWarmupN = Inf; else; endWarmupN = stats.N(idx); end
    
    mu = zeros(1,Niter);
    ln_var = zeros(1,Niter);

    % Parameters for WSABI
    diam = probstruct.PUB - probstruct.PLB;
    kernelCov = diag(diam/10);     % Input length scales for GP likelihood model
    lambda = 1;                     % Ouput length scale for GP likelihood model
    alpha = 0.8;
    
    for iIter = 1:Niter        
        X_train = X(1:N(iIter),:);
        y_train = y(1:N(iIter));
        % Prune trials after warmup
        if N(iIter) >= endWarmupN
            idx_keep = output.X_flag(1:N(iIter));
            X_train = X_train(idx_keep,:);
            y_train = y_train(idx_keep);
        end
        X_iter{iIter} = X_train;
        y_iter{iIter} = y_train;
        lnp = infbench_lnprior(X_train,probstruct);
        y_train = y_train - lnp;  % Remove log prior for WSABI
        [mu(iIter),ln_var(iIter)] = ...
            wsabi_oneshot('L',probstruct.PriorMean,diag(probstruct.PriorVar),kernelCov,lambda,alpha,X_train,y_train);
    end
    vvar = max(real(exp(ln_var)),0);
    
    [history,post] = ...
        StoreAlgoResults(probstruct,[],Niter,X_iter{Niter},y_iter{Niter},mu,vvar,X_iter,y_iter,TotalTime);
    history.scratch.output = output;
    
end

% Remove training data from GPs, too bulky (can be reconstructed)
for i = 1:numel(stats.gp)
    if ~any(stats.funccount(i) == history.SaveTicks)
        stats.gp(i).X = [];
        stats.gp(i).y = [];
    end
end
stats.gpHypFull = [];   % Too bulky
    
history.Output.stats = stats;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gsKL,Mean,Cov,Mode,MTV,xx,Test] = computeStats(vp,gp,probstruct,PosteriorMCMC,VarThresh)
%COMPUTE_STATS Compute additional statistics.
    
% Compute Gaussianized symmetric KL-divergence with ground truth
if PosteriorMCMC == 0
    Ns_moments = 1e6;
    xx = vbmc_rnd(vp,Ns_moments,1,1);
else
    xx = gpsample_vbmc(vp,gp,PosteriorMCMC,1);
end
Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
Nopts = 1 + round(100/vp.K);
Mode = vbmc_mode(vp,Nopts,1);

% Compute marginal total variation
try
    MTV = ComputeMarginalTotalVariation(xx,probstruct);
catch
    MTV = NaN(1,vp.D);
end

% Compute test statistics (problem-specific)
if probstruct.ComputeTestStatistics
    fun = str2func(['infbench_' probstruct.Prob]);
    Test = fun('test',probstruct.ProbInfo,xx);
    Test
else
    Test = [];
end

gsKL
MTV

end
