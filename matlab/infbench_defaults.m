function defaults = infbench_defaults(type,probstruct,options,varargin)
%INFBENCH_DEFAULTS Return default options structure.

if nargin < 2; probstruct = []; end
if nargin < 3; options = []; end

switch lower(type)
    case 'options'
        defaults.OutputDataPrefix = '';  % Prefix to all output files
        defaults.CharFileSep = '@';      % Char separating parts in path and file names

        defaults.PathDirectory = [];     % Path of matlab files
        defaults.RootDirectory = [];     % Root of benchmark results directory tree
        defaults.ProblemDirectory = 'C:/Users/Luigi/Dropbox/Postdoc/BenchMark/problems-data';  % Path of benchmark problem files
        defaults.Display = 'on';         % Level of display
        defaults.ScaleVariables = 'on';  % Center and rescale variables
        defaults.TolFun = 1e-6;          % Required tolerance on function
        defaults.TolX = 1e-6;            % Tolerance on X distance
        defaults.StartFromMode = 0;      % Start from (posterior) mode?
        defaults.OptimizeFirst = 0;      % Run quick optimization first
        defaults.MaxFunEvalMultiplier = 1;      % Increase # func evals
        defaults.StopSuccessfulRuns = 1; % Stop runs when successful
        defaults.SpeedTests = 10;        % Iterations of speed tests
        defaults.ForceFiniteBounds = 0;  % By default do not force finite bounds
        defaults.MaxStoredSamples = 5e3; % Stored posterior MCMC samples
        defaults.Debug = 0;              % Debug mode
    
        defaults.LineStyle = {'-','-','-','-','-.','-.','-.','-.','-.','-.','-.','-','-.','-','-','-','-','-'};
        % defaults.LineStyle = {'-','-.','-','-','-','-','-','-','-','-'};
        defaults.LineColor = [  ...
            141 211 199; ...            % fminsearch
            251 128 114; ...            % fmincon
            128 177 211; ...            % patternsearch
            253 180 98; ...             % mcs
            160 120 100; ...            % global
            70 70 233; ...              % random search
            252 205 229; ...            % simulated annealing
            165 211 195; ...            % genetic algorithm
            159 212 105; ...            % particle swarm
            188 128 189; ...            % cma-es
            212 148 169; ...            % bipop-cma-es            
            88 198 89; ...              % cma-es
            112 248 129; ...            % bipop-cma-es            
            60 220 200; ...            % global
            170 70 133; ...              % random search
            41 111 199; ...            % fminsearch
            151 228 214; ...            % fmincon
            0 0 0 ...                   % bps
            ]/255;
        
    case 'problem'
        
        defaults.MaxFunEvals = min(100 + 50*probstruct.D,750);
        defaults.TolFun = 1e-6;             % Precision data
        defaults.SaveTicks = [];            % Time to save data
        defaults.Noise = [];
        defaults.NoiseSigma = [];           % Added base artificial noise
        defaults.NoiseIncrement = [];       % Added heteroskedastic artificial noise
        defaults.NoiseEstimate = [];        % Estimated size of noise   
        defaults.LocalDataFile = [];        % Local data file to be moved to each local folder
        defaults.VariableComputationTime = false;   % For most problems, computation time is (virtually) constant
        defaults.NonAdmissibleFuncValue = log(realmin);
        defaults.AddLogPrior = false;
        
    case {'plot'}
        defaults.BestOutOf = 1;
        defaults.NumZero = 1e-3;
        defaults.Method = 'FS';
        % defaults.Method = 'IR';
        defaults.SolveThreshold = 10.^(1:-0.1:-2);
%        defaults.VerticalThreshold = 500;
        defaults.Noisy = 0;
        defaults.UnknownMin = 0;
        defaults.DisplayFval = 1;
        defaults.FunEvalsPerD = 500;
        
    case {'plot_noisy'}
        defaults.BestOutOf = 1;
        defaults.NumZero = 1e-2;
        defaults.Method = 'FS';
        defaults.SolveThreshold = 10.^(1:-0.1:-1);
%        defaults.VerticalThreshold = 500;
        defaults.Noisy = 1;
        defaults.UnknownMin = 0;
        defaults.DisplayFval = 1;        
        defaults.FunEvalsPerD = 200;
        
    case 'style'
        name = varargin{1};
        %sep = find(name == '_' | name == '@',1);
        %name(1:sep) = [];
        sep = find(name == '_' | name == '@',1);
        if isempty(sep)
            algo = name;
            algoset = 'base';
        else
            algo = name(1:sep-1);
            algoset = name(sep+1:end);
        end
        line_deterministic = '-';
        line_stochastic = '--';
        
        defaults.color = [0 0 0]/255;
        defaults.linewidth = 2;
        defaults.linestyle = line_deterministic;
        defaults.marker = '';
        defaults.name = [];
        
        switch algo
            case 'vbmc'
                switch algoset
                    % VBMC main paper
                    case {'base','acqproreg'}; defaults.color = [0 0 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = '-'; defaults.name = 'vbmc-P';
                    case {'acqus','acqusreg'}; defaults.color = [180 0 80]/255; defaults.marker = '*'; defaults.linewidth = 3; defaults.linestyle = '-'; defaults.name = 'vbmc-U';
                    case 'acqf2reg'; defaults.color = [0 100 140]/255; defaults.marker = 'o'; defaults.linewidth = 2; defaults.linestyle = '-'; defaults.name = 'vbmc-F';
                    case 'se'; defaults.color = [60 60 60]/255; defaults.marker = '*'; defaults.linewidth = 2; defaults.linestyle = ':';
                    case 'const'; defaults.color = [60 60 60]/255; defaults.marker = 'd'; defaults.linewidth = 2; defaults.linestyle = '--';
                    case 'acqfregvlnn'; defaults.color = [180 0 80]/255; defaults.marker = 'o'; defaults.linewidth = 2; defaults.linestyle = '-.';
                    case {'acqfregvsqrtn','acqmistep1'}; defaults.color = [250 0 80]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = '-';
%                     % VBMC workshop paper                        
%                     case {'base','acqproreg'}; defaults.color = [0 0 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = '-'; defaults.name = 'def';
%                     case {'acqus','acqusreg'}; defaults.color = [180 0 80]/255; defaults.marker = '*'; defaults.linewidth = 2; defaults.linestyle = '-'; defaults.name = 'a-us';
%                     case 'acqf2reg'; defaults.color = [0 100 140]/255; defaults.marker = 'o'; defaults.linewidth = 2; defaults.linestyle = '-'; defaults.name = 'a-gpus';
%                     case 'se'; defaults.color = [60 60 60]/255; defaults.marker = '*'; defaults.linewidth = 2; defaults.linestyle = ':'; defaults.name = 'm-se';
%                     case 'const'; defaults.color = [60 60 60]/255; defaults.marker = 'd'; defaults.linewidth = 2; defaults.linestyle = '--'; defaults.name = 'm-cn';
%                     case 'acqfregvlnn'; defaults.color = [80 0 250]/255; defaults.marker = 'o'; defaults.linewidth = 2; defaults.linestyle = '-'; defaults.name = 'a-ln';
%                     case 'acqfregvsqrtn'; defaults.color = [250 0 80]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = '-'; defaults.name = 'a-sqrt';
                    case 'control'; defaults.color = [120 100 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = '-'; defaults.name = 'vbmc-control';
                    case {'newdef'}; defaults.color = [0 0 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = '-'; defaults.name = 'vbmc-new';
                    case {'acqimiqrnoiserotocorrupthin100','acqimiqrnoiserotonomcmc','acqimiqrnoiseroto','acqimiqrnoiserotoup'}; defaults.color = [180 0 180]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = ':';
                    case {'new2','acqimiqrnoiseup','acqimiqrnoiseup2'}; defaults.color = [180 0 80]/255; defaults.marker = '^'; defaults.linewidth = 2; defaults.linestyle = '-';
                    case {'acqlcb_overhead','lcbnearest_overhead'}; defaults.color = 150*[1 1 1]/255; defaults.marker = ''; defaults.linewidth = 4; defaults.linestyle = '-.';
                    case {'newdef3','oldsettings','gpfast5'}; defaults.color = [120 100 0]/255; defaults.marker = 's'; defaults.linewidth = 2; defaults.linestyle = '-.';
                    case 'test'; defaults.color = [120 100 0]/255; defaults.marker = '>'; defaults.linewidth = 2; defaults.linestyle = ':';
                        % VBMC with noisy functions
                    case {'gpfast6c'}; defaults.color = [180 0 80]/255; defaults.marker = '^'; defaults.linewidth = 2; defaults.linestyle = '-';
                    case {'acqimiqrnoiserotonomcmc2','acqboth'}; defaults.color = [120 100 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = '-';
                    case {'acqmidtstep1_100','ent2','ent2midtstep1mc'}; defaults.color = [250 0 80]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = '-';
                    case {'ent3','ent2midtmcmc50','ent2midtmcmc50step1'}; defaults.color = [250 0 80]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = ':';
                    case {'acqimiqrnoise','acqnoise'}; defaults.color = [180 0 80]/255; defaults.marker = '^'; defaults.linewidth = 2; defaults.linestyle = '--';
                    case {'acqimiqrnoisenegquad','new3','acqimiqrnoiserotocorrup2','roto'}; defaults.color = [120 100 0]/255; defaults.marker = ''; defaults.linewidth = 3; defaults.linestyle = ':';
                        
                    otherwise
                        defaults.color = [0 80 250]/255; defaults.marker = '+'; defaults.linewidth = 2; defaults.linestyle = '-.';
                end
                        
            case 'smc'
                defaults.color = [150 150 150]/255;
                defaults.linewidth = 3;
                defaults.linestyle = ':';
                defaults.marker = '.';                
                
            case 'ais'
                defaults.color = [240 198 114]/255;
                defaults.linewidth = 3;
                defaults.linestyle = ':';
                defaults.marker = 'd';                
                                
            case 'bmc'
                defaults.color = [250 120 150]/255;
                defaults.linewidth = 2;
                defaults.linestyle = '-';
                defaults.marker = 'o';
                 
            case 'wsabi'
                defaults.linestyle = line_deterministic;
                switch algoset
                    case {'base','base2','ip'}; defaults.color = [251 128 114]/255; defaults.marker = 'o'; defaults.name = 'wsabi-L';
                    case {'mm','mm2'}; defaults.color = [114 128 251]/255; defaults.marker = 's'; defaults.name = 'wsabi-M';
                    case {'search'}; defaults.color = [251 128 114]/255; defaults.marker = 'o'; defaults.name = 'wsabi-L'; defaults.linestyle = ':';
                    case {'mmsearch'}; defaults.color = [114 128 251]/255; defaults.marker = 's'; defaults.name = 'wsabi-M'; defaults.linestyle = ':';
                    % case 'actset'; defaults.color = [128 251 114]/255; defaults.marker = 'd';
                end
                defaults.linewidth = 2;

            case 'bbq'
                switch algoset
                    case 'base'; defaults.color = [158 98 159]/255; defaults.marker = 'o'; defaults.name = 'bbq';
                    case 'marginal'; defaults.color = [255 88 169]/255; defaults.marker = '*'; defaults.name = 'bbq*';
                end
                defaults.linewidth = 2;
                defaults.linestyle = line_deterministic;
                
            case {'agp'}
                switch algoset
                    case 'base'; defaults.color = [60 220 200]/255; defaults.marker = 'o'; defaults.name = 'agp';
                    case 'reg'; defaults.color = [255 88 169]/255; defaults.marker = '*'; defaults.name = 'agp (reg)';
                end                 
                defaults.linewidth = 2;
                defaults.linestyle = '--';
                defaults.marker = '+';

            case {'bape'}
                 defaults.color = [128 177 48]/255;
                defaults.linewidth = 2;
                defaults.linestyle = '--';
                defaults.marker = '+';
                defaults.name = 'bape';
                
            case {'parallelgp'}
                switch algoset
                    case 'base'; defaults.color = [60 220 200]/255; defaults.marker = 'o'; defaults.name = 'gp';
                    case 'reg'; defaults.color = [255 88 169]/255; defaults.marker = '*';
                    otherwise; defaults.color = [128 177 48]/255;
                        
                end                 
                                
        end
        
end