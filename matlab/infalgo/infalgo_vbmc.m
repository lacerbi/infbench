function [history,post,algoptions] = infalgo_vbmc(algo,algoset,probstruct)

algoptions = vbmc('all');                   % Get default settings

ControlRunFlag = false;     % Do NOT run in control mode

algoptions.MinFunEvals = probstruct.MaxFunEvals;
algoptions.MaxFunEvals = probstruct.MaxFunEvals;

% VBMC old defaults -- some of these may have changed
algoptions.FunEvalsPerIter = 5;
algoptions.AcqFcn = '@vbmc_acqskl';
algoptions.SearchAcqFcn = '@acqfreg_vbmc';
algoptions.gpMeanFun = 'negquad';
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
algoptions.WarmupNoImproThreshold = Inf; 
algoptions.TolStableWarmup = 15;
algoptions.TolStableExceptions = 1/8;
algoptions.TolStableCount = 40;
algoptions.WarmupCheckMax = false;
algoptions.SGDStepSize = 0.01;
algoptions.RankCriterion = false;
algoptions.EmpiricalGPPrior = true;
algoptions.gpQuadraticMeanBound = false;
algoptions.WarmupOptions = [];
algoptions.WarmupKeepThreshold = '10*nvars';
algoptions.PruningThresholdMultiplier = 1;
algoptions.NSent = @(K) 100*K;
algoptions.NSentFast = @(K) 100*K;
algoptions.NSentFine = @(K) 2^15*K;
algoptions.NSentBoost = [];
algoptions.NSentFastBoost = [];
algoptions.NSentFineBoost = [];
algoptions.ActiveVariationalSamples = 0;
algoptions.FixedMaxMeanGP = false;
algoptions.GPTrainNinit = 1024;
algoptions.GPTrainNinitFinal = 1024;
algoptions.DetEntropyAlpha = 0;
algoptions.ActiveSampleFullUpdate = false;
algoptions.GPTrainInitMethod = 'sobol';
algoptions.VariationalInitRepo = false;
algoptions.MaxIterStochastic = Inf;
algoptions.GPSampleThin = 5;
algoptions.GPTolOpt = 1e-6;
algoptions.GPTolOptMCMC = 0.1;

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

% New VBMC defaults (need to be set manually here)
newdefaults = algoptions;
newdefaults.MinFinalComponents = 50;
newdefaults.WarmupKeepThreshold = '20*nvars';
newdefaults.PruningThresholdMultiplier = @(K) 1/sqrt(K);
newdefaults.gpQuadraticMeanBound = 1;
newdefaults.EmpiricalGPPrior = 0;
newdefaults.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint);
newdefaults.TolStableExcptFrac = 0.2;
newdefaults.TolStableCount = 50;
newdefaults.WarmupCheckMax = true;
newdefaults.SGDStepSize = 0.005;
newdefaults.NSentFine = '@(K) 2^12*K';
newdefaults.NSentFast = 0;
newdefaults.gpMeanFun = 'negquadfix';
newdefaults.GPTrainInitMethod = 'rand';
newdefaults.GPTrainNinitFinal = 64;
newdefaults.MaxIterStochastic = '100*(2+nvars)';
newdefaults.GPSampleThin = 1;


% Options from current problem
switch algoset
    case {0,'debug'}; algoset = 'debug'; algoptions.Debug = true; algoptions.Plot = 'on'; algoptions.FeatureTest = true;
    case {1,'base'}; algoset = 'base';                                                      % Use defaults
    case {2,'acqusreg'}; algoset = 'acqusreg'; algoptions.SearchAcqFcn = @acqusreg_vbmc;    % Vanilla uncertainty sampling
    case {3,'acqproreg'}; algoset = 'acqproreg'; algoptions.SearchAcqFcn = @acqfreg_vbmc;   % Prospective uncertainty sampling
    case {4,'control'}; algoset = 'control'; ControlRunFlag = true;                         % Control experiment
    case {5,'test'}; algoset = 'test'; algoptions.FeatureTest = true;                       % Test feature
    case {6,'narrow'}; algoset = 'narrow'; algoptions.InitDesign = 'narrow';                % Narrow initialization
    case {7,'control2'}; algoset = 'control2'; ControlRunFlag = true;                       % Control experiment, repeated
    case {8,'test2'}; algoset = 'test2'; algoptions.FeatureTest = true;                     % Test feature (second case)

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
    case {51,'acq2'}; algoset = 'acq2'; algoptions.SearchAcqFcn = {@acqmireg_vbmc,@acqfreg_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {52,'acqmif'}; algoset = 'acqmif'; algoptions.SearchAcqFcn = @acqmifreg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {53,'acqmishaping'}; algoset = 'acqmishaping'; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {54,'acq2shaping'}; algoset = 'acq2shaping'; algoptions.SearchAcqFcn = {@acqfreg_vbmc,@acqmireg_vbmc}; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {55,'acq2hedge'}; algoset = 'acq2hedge'; algoptions.AcqHedge = 1; algoptions.SearchAcqFcn = {@acqfreg_vbmc,@acqmireg_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {56,'acq2shapinghedge'}; algoset = 'acq2shapinghedge'; algoptions.AcqHedge = 1; algoptions.SearchAcqFcn = {@acqfreg_vbmc,@acqmireg_vbmc}; algoptions.NoiseShaping = 1; algoptions.NoiseShapingThreshold = '50*nvars'; algoptions.NoiseShapingFactor = 0.2; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {57,'acq2new'}; algoset = 'acq2new'; algoptions.SearchAcqFcn = {@acqmireg_vbmc,@acqfreg_vbmc}; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {60,'step1migps'}; algoset = 'acqstep1migps'; algoptions.FunEvalsPerIter = 1; algoptions.FunEvalStart = 'D'; algoptions.KfunMax = @(N) N; algoptions.SeparateSearchGP = 1; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {61,'lowent'}; algoset = 'lowent'; algoptions.NSentFine = '@(K) 2^12*K'; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {62,'mapgp'}; algoset = 'mapgp'; algoptions = newdefaults; algoptions.FixedMaxMeanGP = true;
    case {63,'mapgpiso'}; algoset = 'mapgpiso'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfixiso';
    case {64,'mapgp2'}; algoset = 'mapgp2'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfix';
    case {65,'mapgp3'}; algoset = 'mapgp3'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainNinit = 256; algoptions.GPTrainInitMethod = 'rand';
    case {66,'mapgp4'}; algoset = 'mapgp4'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand'; algoptions.Plot = 0;
    case {67,'mapgp5'}; algoset = 'mapgp5'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand'; algoptions.Plot = 0; algoptions.SearchAcqFcn = '@acqfreglog_vbmc'; % algoptions.SearchOptimizer = 'bads'; algoptions.SearchMaxFunEvals = '50*nvars';
    case {68,'mapgp2ent1'}; algoset = 'mapgp2ent1'; algoptions = newdefaults; algoptions.gpMeanFun = 'negquadfix'; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;        
    case {69,'gpfast'}; algoset = 'gpfast'; algoptions = newdefaults; algoptions.GPTrainNinitFinal = 64; algoptions.MaxIterStochastic = '100*(2+nvars)';
    case {70,'gpfast2'}; algoset = 'gpfast2'; algoptions = newdefaults; algoptions.GPSampleThin = 1;
    case {71,'gpfast3'}; algoset = 'gpfast3'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-4; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100;
    case {72,'gpfast3up'}; algoset = 'gpfast3up'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-4;
    case {74,'gpfast3b'}; algoset = 'gpfast3b'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-4; algoptions.GPTolOptMCMC = 1e-3; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100;
    case {75,'gpfast5'}; algoset = 'gpfast5'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 2; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-4; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100;
    case {76,'gpfast5up'}; algoset = 'gpfast5up'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-4; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100; algoptions.ActiveSampleFullUpdate = 1;
    case {77,'gpfast5upalpha'}; algoset = 'gpfast5upalpha'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-4; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100; algoptions.ActiveSampleFullUpdate = 1; algoptions.UpdateRandomAlpha = 1;
    case {78,'gpfast3c'}; algoset = 'gpfast3c'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-4; algoptions.GPTolOptMCMC = 1e-3;
    case {79,'gpfast3d'}; algoset = 'gpfast3d'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; 
    case {80,'gpfast3e'}; algoset = 'gpfast3e'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; algoptions.SearchCacheFrac = 0.01; algoptions.StopWarmupReliability = 100;
    case {81,'gpfast4'}; algoset = 'gpfast4'; algoptions = newdefaults; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; algoptions.SearchMaxFunEvals = '200*D'; algoptions.StopWarmupReliability = 100;
    case {82,'gpfast5e'}; algoset = 'gpfast3e'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; algoptions.SearchCacheFrac = 0.01;
        
    % New defaults
    case {100,'newdef'}; algoset = 'newdef'; algoptions = newdefaults;
    case {101,'newdef2'}; algoset = 'newdef2'; algoptions = newdefaults;
    case {102,'newdef3'}; algoset = 'newdef3'; algoptions = newdefaults;

    % Noise
    case {201,'acqf2new'}; algoset = 'acqf2new'; algoptions.Plot = 0; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
        
    % Information-theoretic
    case {301,'oldsettings'}; algoset = 'oldsettings';
%     case {302,'acqmistep1'}; algoset = 'acqmistep1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmireg_vbmc;
%     case {303,'step1'}; algoset = 'step1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1;
%     case {304,'step5'}; algoset = 'step5'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 5;
%     case {305,'acqmistep5'}; algoset = 'acqmistep5'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 5; algoptions.SearchAcqFcn = @acqmireg_vbmc;
%     case {306,'acq2step1'}; algoset = 'acq2step1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = {@acqfreg_vbmc,@acqmireg_vbmc};
%     case {307,'acqmidtstep1_99'}; algoset = 'acqmidtstep1_99'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.RepeatedAcqDiscount = 0.99;
%     case {308,'acqmidtstep1_95'}; algoset = 'acqmidtstep1_99'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.RepeatedAcqDiscount = 0.95;
%     case {309,'acqmidtstep1_100'}; algoset = 'acqmidtstep1_100'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.RepeatedAcqDiscount = 1;
%     case {310,'dtstep1_99'}; algoset = 'dtstep1_99'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqfdtreg_vbmc; algoptions.RepeatedAcqDiscount = 0.99;
%     case {311,'dtstep1_95'}; algoset = 'dtstep1_95'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqfdtreg_vbmc; algoptions.RepeatedAcqDiscount = 0.95;
%     case {312,'dtstep1_100'}; algoset = 'dtstep1_100'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqfdtreg_vbmc; algoptions.RepeatedAcqDiscount = 1;
%     case {321,'step1beta'}; algoset = 'step1beta'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.ELCBOWeight = @(N) sqrt(0.2*2*log(probstruct.D*N^2*pi^2/(6*0.1))); algoptions.Plot = 1;
%     case {322,'step5beta'}; algoset = 'step5beta'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 5; algoptions.ELCBOWeight = @(N) sqrt(0.2*2*log(probstruct.D*N^2*pi^2/(6*0.1))); algoptions.Plot = 1;
%     case {323,'step5K'}; algoset = 'step5K'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 5; algoptions.KfunMax = @(N) min(2,max(2,round(0.5*sqrt(N)))); algoptions.Plot = 1;
%     case {324,'acqmidtstep1K_99'}; algoset = 'acqmidtstep1K_99'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.RepeatedAcqDiscount = 0.99; algoptions.KfunMax = @(N) min(Inf,max(2,floor(0.5*sqrt(N)))); algoptions.Plot = 0;
     case {302,'acqmi'}; algoset = 'acqmi'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc;
     case {303,'acqmiup'}; algoset = 'acqmiup'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1;
     case {304,'acqmiupfast'}; algoset = 'acqmiupfast'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '200*D';
     case {305,'acqmiupfastalpha'}; algoset = 'acqmiupfastalpha'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '200*D'; algoptions.UpdateRandomAlpha = 1;
     case {306,'fast'}; algoset = 'fast'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 2; algoptions.NSgpMaxMain = 2; algoptions.SearchMaxFunEvals = '200*D';
     case {307,'acqmiupfast0'}; algoset = 'acqmiupfast0'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.NSgpMaxWarmup = 2; algoptions.NSgpMaxMain = 2; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; algoptions.StopWarmupReliability = 100; algoptions.SearchCacheFrac = 0.01;
     case {308,'upfast0'}; algoset = 'upfast0'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.SearchMaxFunEvals = '200*D'; algoptions.ActiveSampleFullUpdate = 1; algoptions.GPTolOpt = 1e-4; algoptions.StopWarmupReliability = 100;
     case {309,'acqmiupfast0alpha'}; algoset = 'acqmiupfast0alpha'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.ActiveSampleFullUpdate = 1; algoptions.NSgpMaxWarmup = 0; algoptions.NSgpMaxMain = 0; algoptions.GPTolOpt = 1e-5; algoptions.GPTolOptMCMC = 1e-2; algoptions.StopWarmupReliability = 100; algoptions.SearchCacheFrac = 0.01; algoptions.UpdateRandomAlpha = 1;
                    
    % Entropy tests   
%     case {401,'ent1'}; algoset = 'ent1'; algoptions = newdefaults; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {402,'ent2'}; algoset = 'ent2'; algoptions = newdefaults; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {403,'ent3'}; algoset = 'ent3'; algoptions = newdefaults; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = 0; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {404,'ent2mcmc50'}; algoset = 'ent2mcmc50'; algoptions = newdefaults; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {405,'ent2midtmcmc50'}; algoset = 'ent2midtmcmc50'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {406,'ent2midt'}; algoset = 'ent2midt'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {407,'ent2midtmcmc50step1'}; algoset = 'ent2midtmcmc50step1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {408,'ent2midtstep1'}; algoset = 'ent2midtstep1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {409,'ent1midtstep1'}; algoset = 'ent1midtstep1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {410,'ent2midtmcmc50step1mc'}; algoset = 'ent2midtmcmc50step1mc'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 2; algoptions.NSgpMaxMain = 2; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {411,'ent1midtmcmc50step1'}; algoset = 'ent1midtmcmc50step1'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {412,'ent2midtup'}; algoset = 'ent2midtup'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {413,'ent1midtup'}; algoset = 'ent1midtup'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {414,'ent1midtupmcmc50'}; algoset = 'ent1midtupmcmc50'; algoptions = newdefaults; algoptions.Plot = 1; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveVariationalSamples = 50;
%     case {415,'ent1midtstep1mc'}; algoset = 'ent1midtstep1mc'; algoptions = newdefaults; algoptions.NSgpMaxMain = 2; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {416,'ent2midtstep1mc'}; algoset = 'ent2midtstep1mc'; algoptions = newdefaults; algoptions.NSgpMaxMain = 2; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K;
%     case {417,'ent2midtstep1mapgp'}; algoset = 'ent2midtstep1mapgp'; algoptions = newdefaults; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainNinit = 256; algoptions.GPTrainInitMethod = 'rand';
%     case {418,'ent2mapgp'}; algoset = 'ent2mapgp'; algoptions = newdefaults; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand'; algoptions.Plot = 0;
%     case {419,'ent2midtstep1mapgpbest'}; algoset = 'ent2midtstep1mapgpbest'; algoptions = newdefaults; algoptions.SearchCMAESbest = true; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainNinit = 256; algoptions.GPTrainInitMethod = 'rand';
%     case {420,'ent2midtstep1mapgpfast'}; algoset = 'ent2midtstep1mapgpfast'; algoptions = newdefaults; algoptions.Plot = 0; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '100*D'; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand';
%     case {421,'ent2midtstep1mapgpfastmu'}; algoset = 'ent2midtstep1mapgpfastmu'; algoptions = newdefaults; algoptions.VariableMeans = false; algoptions.Warmup = false; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '200*D'; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSent = 0; algoptions.NSentBoost = @(K) 100*K; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainNinit = 256; algoptions.GPTrainInitMethod = 'rand';
%     case {422,'ent1midtstep1mapgpfast'}; algoset = 'ent1midtstep1mapgpfast'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '100*D'; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand';
     case {423,'up'}; algoset = 'up'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 1;
     case {424,'midtup'}; algoset = 'midtup'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc;
     case {425,'midt'}; algoset = 'midt'; algoptions = newdefaults; algoptions.SearchAcqFcn = @acqmidtreg_vbmc;
%     case {424,'up'}; algoset = 'upalpha0.5'; algoptions = newdefaults; algoptions.Plot = 1; algoptions.gpMeanFun = 'negquadfix'; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.ActiveSampleFullUpdate = 1; algoptions.DetEntropyAlpha = 0.5;
%     case {425,'ent1midtmapgpfastup'}; algoset = 'ent1midtstep1mapgpfastup'; algoptions = newdefaults; algoptions.Plot = 1; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '100*D'; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc;
%     case {431,'ent1midtmapgpfastupa0.5'}; algoset = 'ent1midtstep1mapgpfastup'; algoptions = newdefaults; algoptions.NSgpMaxWarmup = 3; algoptions.NSgpMaxMain = 3; algoptions.SearchMaxFunEvals = '100*D'; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand'; algoptions.DetEntropyAlpha = 0;
%     case {432,'ent1midtmapgpupa0.5'}; algoset = 'ent1midtstep1mapgpfastup'; algoptions = newdefaults; algoptions.ActiveSampleFullUpdate = 1; algoptions.SearchAcqFcn = @acqmidtreg_vbmc; algoptions.NSentFast = 0; algoptions.NSentFastBoost = 0; algoptions.NSentFine = @(K) 2^12*K; algoptions.NSentFineBoost = @(K) 2^12*K; algoptions.gpMeanFun = 'negquadfix'; algoptions.GPTrainInitMethod = 'rand'; algoptions.DetEntropyAlpha = 0;
    
    % Variational active sampling
    case {1000,'vas'}; algoset = 'vas'; 
        
    otherwise
        error(['Unknown algorithm setting ''' algoset ''' for algorithm ''' algo '''.']);
end

% Increase base noise with noisy functions
if ~isempty(probstruct.Noise) || probstruct.IntrinsicNoisy
    algoptions.UncertaintyHandling = 'on';
    NoiseEstimate = probstruct.NoiseEstimate;
    if isempty(NoiseEstimate); NoiseEstimate = 1; end
    algoptions.SpecifyTargetNoise = true;
    % algoptions.NoiseSize = NoiseEstimate(1);
    algoptions.TolStableCount = ceil(algoptions.TolStableCount*1.5);
    algoptions.TolStableWarmup = algoptions.TolStableWarmup*2;
else
    algoptions.UncertaintyHandling = 'off';
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
[vp,elbo,elbo_sd,exitflag,~,output,stats] = ...
    vbmc(@(x) infbench_func(x,probstruct),x0,LB,UB,PLB,PUB,algoptions);
TotalTime = toc(algo_timer);

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
    [post.gsKL,post.Mean,post.Cov,post.Mode] = computeStats(vp,probstruct);

    % Return estimate, SD of the estimate, and gauss-sKL with true moments
    Nticks = numel(history.SaveTicks);
    for iIter = 1:Nticks
        fprintf('%d..',iIter);
        
        idx = find(stats.funccount == history.SaveTicks(iIter),1);
        if isempty(idx); continue; end
        
        % Compute variational solution that VBMC would return at a given iteration
        [vp,~,~,idx_best] = ...
            best_vbmc(stats,idx,algoptions.BestSafeSD,algoptions.BestFracBack,algoptions.RankCriterion);                
        [vp,elbo,elbo_sd] = ...
            finalboost_vbmc(vp,idx_best,[],stats,options_vbmc);
        
        history.Output.N(iIter) = history.SaveTicks(iIter);
        history.Output.lnZs(iIter) = elbo;
        history.Output.lnZs_var(iIter) = elbo_sd^2;
        [gsKL,Mean,Cov,Mode] = computeStats(vp,probstruct);
        history.Output.Mean(iIter,:) = Mean;
        history.Output.Cov(iIter,:,:) = Cov;
        history.Output.gsKL(iIter) = gsKL;
        history.Output.Mode(iIter,:) = Mode;    
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
function [gsKL,Mean,Cov,Mode] = computeStats(vp,probstruct)
%COMPUTE_STATS Compute additional statistics.
    
% Compute Gaussianized symmetric KL-divergence with ground truth
Ns_moments = 1e6;
xx = vbmc_rnd(vp,Ns_moments,1,1);
Mean = mean(xx,1);
Cov = cov(xx);
[kl1,kl2] = mvnkl(Mean,Cov,probstruct.Post.Mean,probstruct.Post.Cov);
gsKL = 0.5*(kl1 + kl2);

% Compute mode
Nopts = 1 + round(100/vp.K);
Mode = vbmc_mode(vp,Nopts,1);

end
