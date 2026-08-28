%==========================================================================
% WANNIERBSE: Bethe-Salpeter Equation Solver for 2D Materials
%
% Version: 1.1.0
% Research Group: Prof. Shun-Jen Cheng (SJCheng Group @ NYCU)
% Technical Contact: https://quantum.web.nycu.edu.tw/wannierbse
%
% DESCRIPTION:
% This is the primary execution script for the WannierBSE package. It
% orchestrates k-mesh generation, tight-binding data ingestion,
% dielectric function processing, and the parallel construction and
% diagonalization of the BSE Hamiltonian (HBSE).
%
% COPYRIGHT c 2026 THE WANNIERBSE DEVELOPMENT TEAM. ALL RIGHTS RESERVED.
%==========================================================================

% --- Initialize Global Environment ---
clear; clc; close all;

t_wbse_total = tic;

%==========================================================================
% [1] SYSTEM DIRECTORY & PATH ARCHITECTURE
%==========================================================================
CurrentPath     = pwd;
ParamPath       = fullfile(CurrentPath, 'Parameters');
UserInputPath   = fullfile(CurrentPath, 'User_input');
CachePath       = fullfile(CurrentPath, 'Precomputed_data');
OutputPath      = fullfile(CurrentPath, 'Exciton_data');
FunctionPath    = fullfile(CurrentPath, 'Functions_for_WBSE');

% Ensure integrity of directory structure
dirs = {ParamPath, UserInputPath, CachePath, OutputPath, FunctionPath, fullfile(CachePath, 'TB_data')};
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir'), mkdir(dirs{i}); end
end

% Register computational kernels
addpath(FunctionPath);
addpath(fullfile(FunctionPath, 'Loaders'));
addpath(fullfile(FunctionPath, 'Utils'));
addpath(fullfile(FunctionPath, 'TB'));
addpath(fullfile(FunctionPath, 'Kmesh'));
addpath(fullfile(FunctionPath, 'Dielectric'));
addpath(fullfile(FunctionPath, 'Wannier_engine'));
addpath(fullfile(FunctionPath, 'Core'));
addpath(fullfile(FunctionPath, 'Core', 'Coulomb'));
addpath(fullfile(FunctionPath, 'Core', 'Direct'));
addpath(fullfile(FunctionPath, 'Core', 'Exchange'));
fprintf('\n[INIT ] WannierBSE Version 1.1.0 environment ready.\n');

%==========================================================================
% [2] PARAMETER INGESTION & GEOMETRY PARSING
%==========================================================================
ReturnPath = pwd;
cd(ParamPath);

% Load simulation control parameters
params = read_params('control.txt');
fn = fieldnames(params);
for k = 1:numel(fn), assignin('base', fn{k}, params.(fn{k})); end

% Load exchange interaction control
exch_ctrl_path = fullfile(ParamPath, 'exchange_control.txt');
if isfile(exch_ctrl_path)
    exch_params = read_params(exch_ctrl_path);
    if isfield(exch_params, 'm_max') && ~isfield(exch_params, 'M_max')
        warning('WBSE:DeprecatedExchangeParameter', ...
                'exchange_control.txt parameter m_max is deprecated; use M_max.');
        exch_params.M_max = exch_params.m_max;
    end
    if isfield(exch_params, 'g_max') && ~isfield(exch_params, 'G_max')
        warning('WBSE:DeprecatedExchangeParameter', ...
                'exchange_control.txt parameter g_max is deprecated; use G_max.');
        exch_params.G_max = exch_params.g_max;
    end
    if isfield(exch_params, 'm_max_tolerance') && ~isfield(exch_params, 'M_max_tolerance')
        warning('WBSE:DeprecatedExchangeParameter', ...
                'exchange_control.txt parameter m_max_tolerance is deprecated; use M_max tolerance.');
        exch_params.M_max_tolerance = exch_params.m_max_tolerance;
    end
    fn_ex = fieldnames(exch_params);
    for k = 1:numel(fn_ex), assignin('base', fn_ex{k}, exch_params.(fn_ex{k})); end

    if ~isfield(exch_params, 'Exchange_Interaction')
        Exchange_Interaction = false;
    end
    if ~isfield(exch_params, 'Use_Precomputed_Chi')
        Use_Precomputed_Chi = false;
    end
else
    Exchange_Interaction = false;
    Use_Precomputed_Chi = false;
end

if ~Exchange_Interaction
    ExchDataPath = fullfile(OutputPath, 'Exchange_data');
    if exist(ExchDataPath, 'dir')
        rmdir(ExchDataPath, 's');
    end
end

% Extract lattice geometry and Wannier Function indices
[a1, a2, a3, a, tau, N_atoms] = read_structure();
cd(ReturnPath);

% User-provided kmesh/TB_data is an external direct-only model path.
% Its TB coefficient columns are assumed to match Parameters/WF_centers.txt.
externalTB = validate_external_kmesh_tb_input(UserInputPath, Nv, Nc, Exchange_Interaction);
useUserProvidedKmeshTB = externalTB.useUserProvidedKmeshTB;
userBandsPath = externalTB.userBandsPath;

% --- Stage 3.4: Wannier Function Preprocessing (Preprocessing Stage) ---
WannierCachePath = fullfile(CachePath, 'Wannier_functions_WBSE');

if useUserProvidedKmeshTB
    % External analytical/TB models keep their own basis convention.
    % Direct interaction only requires coefficient columns to match this file.
    wf_centers_file = fullfile(ParamPath, 'WF_centers.txt');
    fprintf('[LOAD ] Wannier: External TB input detected. Using Parameters/WF_centers.txt without WBSE spin-block reordering.\n');
else
    wf_centers_file = fullfile(CachePath, 'WF_centers_WBSE.txt');
end

% Rule: Volumetric preprocessing is only needed when chi will be calculated.
% Internally generated Wannier90 TB uses the reordered WF_centers_WBSE basis.
if ~useUserProvidedKmeshTB && Exchange_Interaction && ~Use_Precomputed_Chi
    run_wannier_pre = false;
    if ~isfile(wf_centers_file)
        run_wannier_pre = true;
    else
        % Check for the 5 mandatory volumetric cache files
        wannier_files = {'WF_up.mat', 'WF_down.mat', 'r.mat', 'GridInfo.txt', 'WF_AtomPosition.txt'};
        for wf = 1:length(wannier_files)
            if ~isfile(fullfile(WannierCachePath, wannier_files{wf}))
                run_wannier_pre = true; break;
            end
        end
    end

    if run_wannier_pre
        WF_to_mat; % Execute preprocessing module (Granular timers added internally)
    else
        fprintf('[CACHE] Wannier: Found preprocessed data. Skipping...\n');
    end
elseif ~useUserProvidedKmeshTB
    if Exchange_Interaction
        fprintf('[CACHE] Wannier: Precomputed chi requested. Skipping volumetric preprocessing.\n');
    else
        fprintf('[INFO ] Wannier: Exchange Interaction is OFF. Skipping volumetric preprocessing.\n');
    end

    % Internal TB coefficients use the WBSE spin-block basis convention.
    if ~isfile(wf_centers_file)
        fprintf('[GEN  ] Wannier: Generating WF_centers_WBSE.txt...\n');
        wf_centers_in = fullfile(ParamPath, 'WF_centers.txt');
        reorder_WF_centers(wf_centers_in, wf_centers_file);
    end
end

% Load Wannier Function centers and indices
[WF_center, Basis_tau_index] = read_WF_center(N_atoms, wf_centers_file);

% Initialize physical constants and basis sizing
NB = size(WF_center, 1);
[e, eps0] = phys_constant();

%==========================================================================
% [3] DATA ACQUISITION & MULTI-TIER RESOLUTION
%==========================================================================
isExternalKmesh = false; % Internal flag for temporary cache management
kmeshGeneratedThisRun = false;
index_map = [];          % Synchronization vector for Scenario 3
k_fixed = [];            % Sorted physical k-coordinates
kmeshCacheFile = fullfile(CachePath, 'kmesh_WBSE.txt');
tbCachePath = fullfile(CachePath, 'TB_data');
TempTopologyPath = fullfile(CachePath, 'TB_data', 'temp');

% Internal generated k-mesh and TB bands are a coherent cache group.
% Check TB cache completeness before kmesh generation so an incomplete band
% cache does not force a second kmesh pass later in the same run.
userBandsExist = useUserProvidedKmeshTB;
cacheBandsExist = true;
for n = 1:Nv
    if ~isfile(fullfile(tbCachePath, sprintf('v%d_TB_WBSE.txt', n)))
        cacheBandsExist = false; break;
    end
end
if cacheBandsExist
    for n = 1:Nc
        if ~isfile(fullfile(tbCachePath, sprintf('c%d_TB_WBSE.txt', n)))
            cacheBandsExist = false; break;
        end
    end
end

% Resolve kmesh source before dielectric so external user kmesh avoids
% internal cache cleanup, but defer kmesh work until after dielectric.
[kmesh_file, kmesh_tag] = resolve_data_path('kmesh.txt');

% --- Stage 3.1: Dielectric Screening (Screened Coulomb) ---
[eps_file, eps_tag] = resolve_data_path('epsilon.txt');
if strcmp(eps_tag, 'GEN  ')
    fprintf('[GEN  ] Dielectric: Resolving screening function...\n');

    % Load validated SI constant library
    mat_path = fullfile(FunctionPath, 'Utils', 'physical_constant.mat');
    if isfile(mat_path)
        phy_struct = load(mat_path);
    else
        [e_v, eps0_v, hbar_v, m0_v] = phys_constant();
        phy_struct.e = e_v; phy_struct.eps0 = eps0_v;
        phy_struct.hbar = hbar_v; phy_struct.m0 = m0_v;
    end

    % Invoke Semiclassical Poisson Solver
    diel_ctrl = fullfile(ParamPath, 'dielectric_control.txt');
    dielectric_function(diel_ctrl, phy_struct, fullfile(CachePath, 'epsilon_WBSE.txt'));
    eps_file = fullfile(CachePath, 'epsilon_WBSE.txt');
else
    fprintf('[%s] Dielectric: Loaded from %s\n', eps_tag, wbse_display_path(eps_file));
end

% --- Stage 3.2: Momentum Grid (K-Mesh) ---
if ~userBandsExist && ~strcmp(kmesh_tag, 'LOAD ') && ~cacheBandsExist
    fprintf('[CACHE] Incomplete internal kmesh/TB cache. Rebuilding both...\n');
    reset_internal_kmesh_tb_cache(kmeshCacheFile, tbCachePath, TempTopologyPath, true);
    kmesh_file = '';
    kmesh_tag = 'GEN  ';
end

if strcmp(kmesh_tag, 'GEN  ')
    fprintf('[GEN  ] Momentum Grid: Generating sampling mesh...\n');
    cd(CachePath); kmesh_for_BSE(); cd(ReturnPath);
    kmesh_file = kmeshCacheFile;
    kmeshGeneratedThisRun = true;

elseif strcmp(kmesh_tag, 'LOAD ')
    fprintf('[LOAD ] Momentum Grid: Found user input. Validating and Sorting...\n');
    [k_fixed, index_map, h_vals] = kmesh_sorter(kmesh_file);
    isExternalKmesh = true;

    % Cache the validated version
    kmesh_file = kmeshCacheFile;
    fid_fix = fopen(kmesh_file, 'w');
    fprintf(fid_fix, '%d\n%d\n%d\n%d\n%d\n', h_vals);
    for i = 1:size(k_fixed, 1)
        fprintf(fid_fix, '%20.16f    %20.16f     %20.16f\n', k_fixed(i,1), k_fixed(i,2), 0.0);
    end
    fclose(fid_fix);
    kmeshGeneratedThisRun = true;
end

% --- Stage 3.3: Quantum States (Bands & Wavefunctions) ---
% Internal generated k-mesh and TB bands are a coherent cache group.
% If the kmesh was regenerated while a complete TB cache existed, rebuild
% bands so both sides come from the same kmesh_control.
if ~userBandsExist && ~isExternalKmesh
    if kmeshGeneratedThisRun && cacheBandsExist
        fprintf('[CACHE] Bands: Existing TB_data does not match regenerated kmesh. Rebuilding...\n');
        reset_internal_kmesh_tb_cache(kmeshCacheFile, tbCachePath, TempTopologyPath, false);
        cacheBandsExist = false;
    end
end

if userBandsExist
    fprintf('[LOAD ] Bands: Found user input in %s. Loading and synchronizing...\n', wbse_display_path(userBandsPath));
    if isExternalKmesh
        DestTB = fullfile(CachePath, 'TB_data');
        if ~exist(DestTB, 'dir'), mkdir(DestTB); end
        kmesh_band_sync(userBandsPath, DestTB, index_map, Nv, Nc, NB, k_fixed);
        bandSource = DestTB;
    else
        bandSource = userBandsPath;
    end
elseif cacheBandsExist
    fprintf('[LOAD ] Bands: Found precomputed bands. Loading...\n');
    bandSource = tbCachePath;
else
    fprintf('[GEN  ] Bands: Resolving Wannier Tight-Binding Hamiltonian...\n');
    tbContext = struct();
    tbContext.kmeshFile = kmesh_file;
    tbContext.wannierWinFile = fullfile(UserInputPath, 'wannier90.win');
    tbContext.structureFile = fullfile(UserInputPath, 'structure.txt');
    tbContext.wannierHrFile = fullfile(UserInputPath, 'wannier90_hr.dat');
    tbContext.wfCentersFile = fullfile(ParamPath, 'WF_centers.txt');
    tbContext.controlFile = fullfile(ParamPath, 'control.txt');
    tbContext.wtbControlFile = fullfile(ParamPath, 'WTB_control.txt');
    tbContext.tbOutputDir = tbCachePath;
    tbContext.exchangeInteraction = Exchange_Interaction;
    TB_for_BSE(tbContext);
    bandSource = tbCachePath;
end

% Topology sidecars are only required while generating TB bands.
if isfolder(TempTopologyPath)
    rmdir(TempTopologyPath, 's');
end

%==========================================================================
% [4] NUMERICAL SOLVER INITIALIZATION
%==========================================================================
[q_data, epsQ2D_data] = read_dielectric(eps_file);
[f, Nv_coarse, Nv_dense, N_co, N_de, k] = read_kmesh(kmesh_file);
N_kpt = size(k, 1);
[b1, b2] = reciprocal_lattice_vector(a1, a2);

% --- High-Performance Computing (HPC) Setup ---
exchange_header_printed = false;
poolobj = gcp('nocreate');
if isempty(poolobj) || poolobj.NumWorkers ~= Npar_HBSE
    if Exchange_Interaction && ~Use_Precomputed_Chi
        fprintf('[PRE  ] Exchange Interaction: Preparing short-range electron-hole Coulomb exchange...\n');
        fprintf('  -> Starting parallel pool for exchange and BSE solver stages: %d workers\n', Npar_HBSE);
        exchange_header_printed = true;
    else
        fprintf('  -> Starting parallel pool for BSE solver kernels: %d workers\n', Npar_HBSE);
    end
    if ~isempty(poolobj), delete(poolobj); end
    parpool('local', Npar_HBSE);
end

% Ingest Quantum Coefficients
[Ev, Cv] = read_VB(bandSource, FunctionPath, Nv, N_kpt, e, NB);
[Ec, Cc] = read_CB(bandSource, FunctionPath, Nc, N_kpt, e, NB);

%==========================================================================
% [5] BSE HAMILTONIAN (HBSE) CONSTRUCTION
%==========================================================================
[Wq0, wi] = singularity_and_kweight(N_de, N_co, Nv_coarse, Nv_dense, d, q_data, epsQ2D_data, a, f);

% --- Short-Range Exchange Pre-calculation ---
Time_rho = NaN;
Time_exchange_tensor = NaN;
if Exchange_Interaction
    t_exchange_pre = tic;
    PrecomputedExchPath = fullfile(CachePath, 'Exchange_interaction');
    chi_cache_file = fullfile(PrecomputedExchPath, 'chi_WBSE.mat');
    R_mNN_cache_file = fullfile(PrecomputedExchPath, 'R_mNN_WBSE.mat');

    if Use_Precomputed_Chi
        fprintf('[LOAD ] Exchange Interaction: Loading precomputed chi/R_mNN tensors...\n');
        fprintf('  -> Source: %s\n', wbse_display_path(PrecomputedExchPath));
        if ~isfile(chi_cache_file)
            error('WBSE:MissingPrecomputedChi', 'Precomputed chi file not found: %s', chi_cache_file);
        end
        if ~isfile(R_mNN_cache_file)
            error('WBSE:MissingPrecomputedRmNN', 'Precomputed R_mNN file not found: %s', R_mNN_cache_file);
        end

        chi_data = load(chi_cache_file, 'chi');
        R_mNN_data = load(R_mNN_cache_file, 'R_mNN');
        if ~isfield(chi_data, 'chi')
            error('WBSE:InvalidPrecomputedChi', 'File does not contain variable chi: %s', chi_cache_file);
        end
        if ~isfield(R_mNN_data, 'R_mNN')
            error('WBSE:InvalidPrecomputedRmNN', 'File does not contain variable R_mNN: %s', R_mNN_cache_file);
        end

        chi = chi_data.chi;
        R_mNN = R_mNN_data.R_mNN;
        if size(chi, 5) ~= size(R_mNN, 1) || size(chi, 6) ~= size(R_mNN, 1)
            error('WBSE:IncompatiblePrecomputedExchange', ...
                  'Precomputed chi and R_mNN sizes are incompatible.');
        end

        chi_p = permute(chi, [3, 4, 1, 2, 5, 6]);
        fprintf('  -> Precomputed exchange tensors loaded and validated.\n');
        clear chi chi_data R_mNN_data;
    else
        if ~exchange_header_printed
            fprintf('[PRE  ] Exchange Interaction: Preparing short-range electron-hole Coulomb exchange...\n');
        end
        fprintf('  -> Loading Wannier functions and real-space grid...\n');

        t_exch_data = tic;
        load(fullfile(WannierCachePath, 'WF_up.mat'));
        load(fullfile(WannierCachePath, 'WF_down.mat'));
        load(fullfile(WannierCachePath, 'r.mat'));

        [Na1_full, Na2_full, Na3_full, deltaV, supercell_size, shift_a1, shift_a2, shift_a3] = read_grid_info(WannierCachePath, a1);
        [J, hu, hv, hz, u, v, u_point, v_point, z_point] = CoordinateTransformation(r, a1, a2, Na3_full);

        [R_mNN, n1_n2_mNN, d_max] = find_rgrid_in_mmn(a1, a2, tau, N_atoms, supercell_size, M_max, nn_distance_tolerance, M_max_tolerance);
        [W_up_t, W_down_t, tin_all, r_t, u_t, v_t, z_t, u_point_t, v_point_t, z_point_t, Nu_t, Nv_t, Nz_t, Na1_t, Na2_t, zcut_min, zcut_max] = ...
                                                       find_truncated_wfs(W_up, W_down, r, u, v, deltaV, Na1_full, z_point, tol_probability);
        fprintf('  -> Real-space grid and Wannier cutoff region prepared in %.2f s\n', toc(t_exch_data));

        % Momentum Mappings (Hardcoded kex=0)
        kex = [0, 0, 0];
        [G_all] = generate_g_vectors(G_max, b1, b2);

        % Integration Kernels (Parallelized)
        fprintf('  -> Calculating rho overlap density...\n');
        t_rho = tic;
        [rho] = rho_calculation(W_up, W_down, W_up_t, W_down_t, G_all, R_mNN, n1_n2_mNN, ...
                                zcut_min, zcut_max, Na1_full, Na2_full, Na3_full, shift_a1, shift_a2, shift_a3, tau, ...
                                N_atoms, d_max, Basis_tau_index, M_max_tolerance, tin_all, Nz_t, Na1_t, Na2_t, ...
                                kex, a1, a2, u_point_t, v_point_t, u_t, hu, hv, J, 0.5, Npar_rho);
        Time_rho = toc(t_rho);
        fprintf('  -> rho overlap density completed in %.2f s\n', Time_rho);

        fprintf('  -> Calculating chi exchange tensor...\n');
        t_chi = tic;
        [chi] = chi_calculation(rho, R_mNN, G_all, kex, tau, N_atoms, d_max, Basis_tau_index, M_max_tolerance, Npar_Xi, e, eps0, z_point_t, b1);
        Time_exchange_tensor = toc(t_chi);
        fprintf('  -> chi exchange tensor completed in %.2f s\n', Time_exchange_tensor);
        chi_p = permute(chi, [3, 4, 1, 2, 5, 6]);

        % Save intermediate exchange interaction data
        ExchDataPath = fullfile(OutputPath, 'Exchange_data');
        if ~exist(ExchDataPath, 'dir'), mkdir(ExchDataPath); end
        if ~exist(PrecomputedExchPath, 'dir'), mkdir(PrecomputedExchPath); end
        save(fullfile(ExchDataPath, 'chi.mat'), 'chi', '-v7.3');
        save(fullfile(ExchDataPath, 'rho.mat'), 'rho', '-v7.3');
        save(fullfile(ExchDataPath, 'G_all.mat'), 'G_all', '-v7.3');
        save(fullfile(ExchDataPath, 'R_mNN.mat'), 'R_mNN', '-v7.3');
        save(chi_cache_file, 'chi', '-v7.3');
        save(R_mNN_cache_file, 'R_mNN', '-v7.3');

        clear W_up W_down rho chi; % Memory optimization
        fprintf('[DONE ] Exchange Interaction: Short-range exchange preprocessing completed in %.2f s\n', toc(t_exchange_pre));
    end
else
    chi_p = 0; R_mNN = 0;
end

vc_band_BSE = generate_vc_band(Nv, Nc);
HBSE = zeros(Nc * Nv * N_kpt, Nc * Nv * N_kpt);

fprintf('[SOLVE] Constructing HBSE matrix...\n');
tic;

% Iterative block-matrix assembly
for vvcci = 1:size(vc_band_BSE, 1)
    vv = vc_band_BSE(vvcci, 1); cc = vc_band_BSE(vvcci, 2);
    fprintf('  -> Processing subspace: v%d, c%d\n', vv, cc);

    for vci = vvcci:size(vc_band_BSE, 1)
        v = vc_band_BSE(vci, 1); c = vc_band_BSE(vci, 2);
        HBSE_element = zeros(N_kpt, N_kpt);

        parfor (j = 1:N_kpt, Npar_HBSE)
            % Kinetic Energy Diagonal contribution
            KE_vec = zeros(N_kpt, 1);
            val_KE = (Ec(j, c) - Ev(j, v)) * (eq(c, cc) * eq(v, vv));
            KE_vec(j) = val_KE;

            % Momentum-dependent Coulomb Kernel
            q = bsxfun(@minus, k, k(j, :));
            [Vq, epsilon_q] = screened_coulomb_interaction(q, b1, b2, d, q_data, epsQ2D_data);

            % Direct Interaction integration
            Vd = direct_coulomb_interaction(Vq, epsilon_q, k, b1, b2, Wq0(j), wi(j), q, tau, Basis_tau_index, WF_center, N_atoms, Cc, Cv, c, cc, v, vv, j, e, eps0);

            % Exchange Interaction (Short-Range)
            Vx = zeros(N_kpt, 1);
            if Exchange_Interaction
                Vx = exchange_coulomb_interaction(chi_p, R_mNN, k, N_kpt, Cv, Cc, c, cc, v, vv, j);
            end

            HBSE_element(:, j) = KE_vec + (wi(j) / (2 * pi)^2) * (Vd + Vx);
        end

        % Map local blocks to global HBSE
        idx1 = (1:N_kpt) + (((c - 1) + (v - 1) * Nc) * N_kpt);
        idx2 = (1:N_kpt) + (((cc - 1) + (vv - 1) * Nc) * N_kpt);
        HBSE(idx1, idx2) = HBSE_element;

        % Symmetrize off-diagonal blocks
       if v ~= vv || c ~= cc
            parfor (j = 1:N_kpt, Npar_HBSE)
                HBSE_element(:, j) = HBSE_element(:, j) / (wi(j) / (2 * pi)^2);
            end
            HBSE(idx2, idx1) = (transpose(wi) / (2 * pi)^2) .* ctranspose(HBSE_element);
        end
    end
end
Time_HBSE = toc;
fprintf('[SOLVE] Computation time of constructing BSE matrix: %.2f s\n', Time_HBSE);


%==========================================================================
% [6] EIGENSOLVER
%==========================================================================
fprintf('[SOLVE] Diagonalizing HBSE matrix...\n');
tic;
[A, E] = eigs(HBSE, NVAL, 'sr', 'Subspacedimension', 2 * NVAL + 20);
Ex = real(diag(E));
Time_eigs = toc;
fprintf('[SOLVE] Computation time of diagonalizing BSE matrix: %.2f s\n', Time_eigs);

% Post-processing of wavefunctions
A = A_output(A, Nc, Nv, wi, N_kpt);

%==========================================================================
% [7] Exciton Binding Energy
%==========================================================================
tic;
vertical_gaps = Ec(:, 1) - Ev(:, 1);
[min_Eg, min_idx] = min(vertical_gaps);
k_at_gap = k(min_idx, :);
gap_time = toc;
% fprintf('[INFO ] Found Direct Gap: %.4f eV at k = [%.3f, %.3f] (Time: %.4f s)\n', min_Eg, k_at_gap(1), k_at_gap(2), gap_time);

%==========================================================================
% [8] DATA EXPORT & SYSTEM CLEANUP
%==========================================================================
Total_Time = Time_HBSE + Time_eigs;
fprintf('[DONE ] Total BSE stage completed in %.2f s\n', Total_Time);

% Save in binary .mat format only
save(fullfile(OutputPath, 'A.mat'), 'A', '-v7.3');
save(fullfile(OutputPath, 'Ex.mat'), 'Ex', '-v7.3');
save(fullfile(OutputPath, 'wi.mat'), 'wi', '-v7.3');

k_save = k; k_save(:, 3) = 0;
save(fullfile(OutputPath, 'k.mat'), 'k_save', '-v7.3');

log_info = struct();
log_info.N_kpt = N_kpt;
log_info.Nv = Nv;
log_info.Nc = Nc;
log_info.Exchange_Interaction = Exchange_Interaction;
log_info.Use_Precomputed_Chi = Use_Precomputed_Chi;
log_info.Npar_TB = Npar_TB;
log_info.Npar_rho = Npar_rho;
log_info.Npar_Xi = Npar_Xi;
log_info.Npar_HBSE = Npar_HBSE;
log_info.Time_rho = Time_rho;
log_info.Time_exchange_tensor = Time_exchange_tensor;
log_info.Time_HBSE = Time_HBSE;
log_info.Time_eigs = Time_eigs;
log_info.Time_BSE_total = Total_Time;
if Exchange_Interaction
    log_info.M_max = M_max;
    log_info.G_max = G_max;
    log_info.tol_probability = tol_probability;
    log_info.M_max_tolerance = M_max_tolerance;
    log_info.nn_distance_tolerance = nn_distance_tolerance;
end

if isExternalKmesh
    if userBandsExist
        ArchiveBands = fullfile(OutputPath, 'TB_Bands');
        if ~exist(ArchiveBands, 'dir'), mkdir(ArchiveBands); end
        movefile(fullfile(CachePath, 'TB_data', 'v*_TB_WBSE.txt'), ArchiveBands);
        movefile(fullfile(CachePath, 'TB_data', 'c*_TB_WBSE.txt'), ArchiveBands);
        archived_files = dir(fullfile(ArchiveBands, '*_WBSE.txt'));
        for f = 1:length(archived_files)
            old_name = archived_files(f).name;
            new_name = strrep(old_name, '_WBSE', '');
            movefile(fullfile(ArchiveBands, old_name), fullfile(ArchiveBands, new_name));
        end
    end
    temp_cache = fullfile(CachePath, 'kmesh_WBSE.txt');
    if isfile(temp_cache), delete(temp_cache); end
end

delete(gcp('nocreate'));
log_info.Time_WBSE_total = toc(t_wbse_total);
output_log(OutputPath, log_info);
fprintf('[DONE ] Total WBSE execution time: %.2f s\n', log_info.Time_WBSE_total);
fprintf('[DONE ] Execution Successful. Output saved to %s/.\n', wbse_display_path(OutputPath));

function displayPath = wbse_display_path(pathValue)
anchors = {'Precomputed_data', 'User_input', 'Parameters', 'Exciton_data'};
displayPath = char(pathValue);
normalizedPath = strrep(displayPath, '\', '/');
for iAnchor = 1:numel(anchors)
    anchor = anchors{iAnchor};
    matchIndex = strfind(normalizedPath, anchor);
    if ~isempty(matchIndex)
        displayPath = normalizedPath(matchIndex(1):end);
        return
    end
end
displayPath = normalizedPath;
end
