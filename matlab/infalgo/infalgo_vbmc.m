function [history,post,algoptions] = infalgo_vbmc(algo,algoset,probstruct)

algoptions = vbmc('all');                   % Get default settings

ControlRunFlag = false;     % Do NOT run in control mode

algoptions.MinFunEvals = probstruct.MaxFunEvals;
algoptions.MaxFunEvals = probstruct.MaxFunEvals;

% VBMC old defaults -- some of these may have changed
algoptions.FunEvalsPerIter = 5;
algoptions.AcqFcn = '@vbmc_acqskl';
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
algoptions.TolStableExceptions = 1/8;
algoptions.TolStableCount = 40;
algoptions.WarmupCheckMax = false;
algoptions.SGDStepSize = 0.01;
algoptions.RankCriterion = false;
algoptions.EmpiricalGPPrior = true;
algoptions.gpQuadraticMeanBound = false;
algoptions.WarmupOptions = [];
algoptions.WarmupKeepThreshold = '10*nvars';

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
    case {58,'step1mi'}; algoset = 'acqstep1mi'; algoptions.FunEvalsPerIter = 1; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {59,'step1'}; algoset = 'acqstep1'; algoptions.FunEvalsPerIter = 1; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {60,'step1migps'}; algoset = 'acqstep1migps'; algoptions.FunEvalsPerIter = 1; algoptions.FunEvalStart = 'D'; algoptions.KfunMax = @(N) N; algoptions.SeparateSearchGP = 1; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
    case {61,'lowent'}; algoset = 'lowent'; algoptions.NSentFine = '@(K) 2^12*K'; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
        
        
    % New defaults
    case {100,'newdef'}; algoset = 'newdef'; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {101,'newdef2'}; algoset = 'newdef2'; algoptions.Plot = 1; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {102,'newdef3'}; algoset = 'newdef3'; algoptions.Plot = 0; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;

    % Noise
    case {201,'acqf2new'}; algoset = 'acqf2new'; algoptions.Plot = 0; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;        
        
    % Information-theoretic
    case {301,'oldsettings'}; algoset = 'oldsettings';
    case {302,'finalK'}; algoset = 'finalK'; algoptions.MinFinalComponents = 60; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {303,'finalKmi'}; algoset = 'finalK'; algoptions.SearchAcqFcn = @acqmireg_vbmc; algoptions.MinFinalComponents = 60; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {304,'prune20'}; algoset = 'prune20'; algoptions.WarmupKeepThreshold = '20*nvars'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
    case {305,'prune50'}; algoset = 'prune50'; algoptions.WarmupKeepThreshold = '50*nvars'; algoptions.gpQuadraticMeanBound = 1; algoptions.EmpiricalGPPrior = 0; algoptions.WarmupNoImproThreshold = 20 + 5*numel(probstruct.InitPoint); algoptions.TolStableExcptFrac = 0.2; algoptions.TolStableCount = 50; algoptions.WarmupCheckMax = true; algoptions.SGDStepSize = 0.005;
        
        
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
    algoptions.TolStableIters = ceil(algoptions.TolStableIters*1.5);
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
