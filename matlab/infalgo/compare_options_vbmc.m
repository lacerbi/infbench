%COMPARE_OPTIONS_VBMC
%
% Compare the current VBMC option defaults with the infbench

if ~exist('algoset','var') || isempty(algoset)
    algoset = 'base';
end

% Add utility folders to path
addpath([fileparts(which('vbmc.m')),filesep,'misc']);
addpath([fileparts(which('vbmc.m')),filesep,'utils']);
addpath([fileparts(which('infbench_run.m')),filesep,'infalgo']);

% Get current infbench options
[options_infbench,probstruct] = infalgo_vbmc('vbmc',algoset,[],true);
nvars = probstruct.D;
options_infbench = setupoptions_vbmc(nvars,options_infbench,options_infbench);

% Get VBMC defaults
options_vbmc = vbmc('all');
options_vbmc = setupoptions_vbmc(nvars,options_vbmc,options_vbmc);

dif_count = 0;

% Go over all fields in the VBMC options and print differences
for ff = fields(options_vbmc)'

    opt_vbmc = options_vbmc.(ff{:});
    opt_infbench = options_infbench.(ff{:});

    if isa(opt_vbmc,'function_handle')
        opt_vbmc = func2str(opt_vbmc);
    end
    if isa(opt_infbench,'function_handle')
        opt_infbench = func2str(opt_infbench);
    end

    try
        if isstring(opt_vbmc) || isstring(opt_infbench)
            isdif = ~strcmp(opt_vbmc, opt_infbench);
        else
            isdif = opt_vbmc ~= opt_infbench;
        end
    catch
        isdif = true;
    end

    if isdif
        dif_count = dif_count + 1;

        fprintf('%s:', ff{:})
        if ~isstruct(opt_vbmc)
            opt_vbmc
            opt_infbench
        end
        fprintf('\n######\n');
    end
end

fprintf('SUMMARY: %d different fields out of %d.\n', dif_count, numel(fields(options_vbmc)));