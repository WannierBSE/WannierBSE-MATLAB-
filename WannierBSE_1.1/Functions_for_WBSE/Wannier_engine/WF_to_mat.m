%==========================================================================
% WF_TO_MAT
%
% DESCRIPTION:
%   Orchestrates Wannier-function preprocessing by importing spinor XSF data, assembling complex components, shifting grids, normalizing wavefunctions, and exporting WBSE-ready caches.
%==========================================================================

fprintf('[PRE  ] Wannier: Preparing real-space Wannier functions...\n');
t_wan_total = tic;

% 1. Initialize Paths
Paths = init_wannier_paths();

% 2. Read User Parameters from exchange_control.txt
if evalin('base', 'exist(''ParamPath'', ''var'')')
    ParamPath = evalin('base', 'ParamPath');
    param_file = fullfile(ParamPath, 'exchange_control.txt');
else
    param_file = fullfile(Paths.parameters, 'exchange_control.txt');
end
user_params = read_params(param_file);

% Extract the cleanup flag (default to true if not found)
if isfield(user_params, 'cleanup_xsf')
    cleanup_xsf = user_params.cleanup_xsf;
else
    cleanup_xsf = true; 
end

% 4. Create the Workspace
if ~exist(Paths.workspace, 'dir')
    mkdir(Paths.workspace);
end

% 5. Execution Pipeline
% (The scripts below will now have access to 'Paths' and 'cleanup_xsf')
fprintf('  -> Importing XSF volumetric data...\n');
t_step = tic;
run('test_up_real.m');
run('test_up_imag.m');
run('test_down_real.m');
run('test_down_imag.m');
fprintf('  -> XSF volumetric data imported in %.2f s\n', toc(t_step));

run('Delete_xsf.m');

t_step = tic;
run('Step2_SaveWF_spinor_v2.m');
fprintf('  -> Spinor basis assembled in %.2f s\n', toc(t_step));

t_step = tic;
run('Step3_shift_r_and_z.m');
fprintf('  -> Spatial grid shifted in %.2f s\n', toc(t_step));

t_step = tic;
run('Step4_Normalize_WFs.m');
fprintf('  -> Wavefunctions normalized in %.2f s\n', toc(t_step));

t_step = tic;
run('Step5_reorder_WFs_TMD_ML_only.m');
fprintf('  -> Wannier basis reordered for the active material convention in %.2f s\n', toc(t_step));

fprintf('[DONE ] Wannier: Preprocessing completed in %.2f s. Cleanup flag: %d\n', toc(t_wan_total), cleanup_xsf);
