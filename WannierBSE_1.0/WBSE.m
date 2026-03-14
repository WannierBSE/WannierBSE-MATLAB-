%==========================================================================
% WANNIERBSE: Bethe-Salpeter Equation Solver for 2D Materials
% 
% Version: 1.0.0 
% Research Group: Prof. Shun-Jen Cheng (SJCheng Group @ NYCU)
% Technical Contact: https://quantum.web.nycu.edu.tw/wannierbse
% 
% DESCRIPTION:
% This is the primary execution script for the WannierBSE package. It 
% orchestrates k-mesh generation, tight-binding data ingestion, 
% dielectric function processing, and the parallel construction and 
% diagonalization of the BSE Hamiltonian (HBSE).
%
% COPYRIGHT © 2026 THE WANNIERBSE DEVELOPMENT TEAM. ALL RIGHTS RESERVED.
%==========================================================================

% --- Initialize Global Environment ---
clear; clc; close all;

%==========================================================================
% [1] SYSTEM DIRECTORY & PATH ARCHITECTURE
%==========================================================================
CurrentPath     = pwd;
ParamPath       = fullfile(CurrentPath, 'Parameters');
UserInputPath   = fullfile(CurrentPath, 'User_input');
CachePath       = fullfile(CurrentPath, 'Precomputed_data');
OutputPath      = fullfile(CurrentPath, 'Exciton_data');
FunctionPath    = fullfile(CurrentPath, 'Functions_for_BSE');

% Ensure integrity of directory structure
dirs = {ParamPath, UserInputPath, CachePath, OutputPath, FunctionPath};
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir'), mkdir(dirs{i}); end
end

% Register computational kernels
addpath(FunctionPath);
fprintf('[INIT ] WannierBSE Version 1.0.0 environment ready.\n');

%==========================================================================
% [2] PARAMETER INGESTION & GEOMETRY PARSING
%==========================================================================
ReturnPath = pwd; 
cd(ParamPath);

% Load simulation control parameters
params = read_params('control.txt'); 
fn = fieldnames(params);
for k = 1:numel(fn), assignin('base', fn{k}, params.(fn{k})); end

% Extract lattice geometry and Wannier Function indices
[a1, a2, a3, a, tau, N_atoms] = Read_win(); 
[WF_center, Basis_tau_index] = Read_WF_center(N_atoms);
cd(ReturnPath);

% Initialize physical constants and basis sizing
NB = size(WF_center, 1); 
[e, eps0] = phys_constant();

%==========================================================================
% [3] DATA ACQUISITION & MULTI-TIER RESOLUTION
%==========================================================================
isExternalKmesh = false; % Internal flag for temporary cache management
index_map = [];          % Synchronization vector for Scenario 3
k_fixed = [];            % Sorted physical k-coordinates

% --- Stage 3.1: Momentum Grid (K-Mesh) ---
[kmesh_file, kmesh_tag] = resolve_data_path('kmesh.txt');

if strcmp(kmesh_tag, 'GEN  ')
    fprintf('[GEN  ] Momentum Grid: Generating sampling mesh...\n');
    cd(CachePath); kmesh_for_BSE(); cd(ReturnPath);
    kmesh_file = fullfile(CachePath, 'kmesh_WBSE.txt');
    
elseif strcmp(kmesh_tag, 'LOAD ')
    fprintf('[LOAD ] Momentum Grid: Found user input. Validating and Sorting...\n');
    [k_fixed, index_map, h_vals] = kmesh_sorter(kmesh_file);
    isExternalKmesh = true; 
    
    % Cache the validated version
    kmesh_file = fullfile(CachePath, 'kmesh_WBSE.txt');
    fid_fix = fopen(kmesh_file, 'w');
    fprintf(fid_fix, '%d\n%d\n%d\n%d\n%d\n', h_vals);
    for i = 1:size(k_fixed, 1)
        fprintf(fid_fix, '%20.16f    %20.16f     %20.16f\n', k_fixed(i,1), k_fixed(i,2), 0.0);
    end
    fclose(fid_fix);
end

% --- Stage 3.2: Dielectric Screening (Screened Coulomb) ---
[eps_file, eps_tag] = resolve_data_path('epsilon.txt');
if strcmp(eps_tag, 'GEN  ')
    fprintf('[GEN  ] Dielectric: Resolving screening function...\n');
    
    % Load validated SI constant library
    mat_path = fullfile(FunctionPath, 'physical_constant.mat');
    if isfile(mat_path)
        phy_struct = load(mat_path);
    else
        [e_v, eps0_v, hbar_v, m0_v] = phys_constant();
        phy_struct.e = e_v; phy_struct.eps0 = eps0_v; 
        phy_struct.hbar = hbar_v; phy_struct.m0 = m0_v;ding Hamiltonians
    end
    
    % Invoke Semiclassical Poisson Solver
    diel_ctrl = fullfile(ParamPath, 'dielectric_control.txt');
    dielectric_function(diel_ctrl, phy_struct, fullfile(CachePath, 'epsilon_WBSE.txt'));
    eps_file = fullfile(CachePath, 'epsilon_WBSE.txt');
else
    fprintf('[%s] Dielectric: Loaded from %s\n', eps_tag, eps_file);
end

% --- Stage 3.3: Quantum States (Bands & Wavefunctions) ---
% Check for existence of external user-provided electronic structure files
userBandsExist = (exist(fullfile(UserInputPath, 'v1_TB.txt'), 'file') || ...
                  exist(fullfile(UserInputPath, 'c1_TB.txt'), 'file'));

% Check for existence of validated precomputed data in cache
cacheBandsExist = (exist(fullfile(CachePath, 'v1_TB_WBSE.txt'), 'file') || ...
                   exist(fullfile(CachePath, 'c1_TB_WBSE.txt'), 'file'));

if userBandsExist
    % ---------------------------------------------------------------------
    % CASE A: EXTERNAL USER DATA (Requires Synchronization)
    % ---------------------------------------------------------------------
    fprintf('[LOAD ] Bands: Found user input. Loading and synchronizing...\n');
    
    if isExternalKmesh
        % Perform block-wise realignment with the sorted momentum grid
        kmesh_band_sync(UserInputPath, CachePath, index_map, Nv, Nc, NB, k_fixed);
        bandSource = CachePath;
    else
        % No k-mesh sorting required; ingest directly from user folder
        bandSource = UserInputPath;
    end

elseif cacheBandsExist
    % ---------------------------------------------------------------------
    % CASE B: PRECOMPUTED CACHE (Validated session data)
    % ---------------------------------------------------------------------
    fprintf('[LOAD ] Bands: Found precomputed bands. Loading...\n');
    bandSource = CachePath;

else
    % ---------------------------------------------------------------------
    % CASE C: INTERNAL GENERATION (Wannier Tight-Binding Solver)
    % ---------------------------------------------------------------------
    fprintf('[GEN  ] Bands: Resolving Wannier Tight-Binding Hamiltonian...\n');
    fprintf('-------------------------------------------------------\n');        
    
    % Initialize transient input environment
    TempInput = fullfile(CurrentPath, 'input_files'); 
    if ~exist(TempInput,'dir'), mkdir(TempInput); end
    
    % Stage necessary Hamiltonian data for internal TB resolution
    if exist(fullfile(UserInputPath, 'wannier90_hr.dat'), 'file')
        copyfile(fullfile(UserInputPath, 'wannier90_hr.dat'), ...
                 fullfile(TempInput, 'wannier90_hr.dat'));
    end
    
    % Execution of the WannierBSE Tight-Binding Engine
    cd(CachePath); 
    try
        TB_for_BSE(); 
    catch ME
        cd(ReturnPath); rmdir(TempInput,'s'); rethrow(ME); 
    end
    cd(ReturnPath); rmdir(TempInput, 's');
    
    bandSource = CachePath;
end

%==========================================================================
% [4] NUMERICAL SOLVER INITIALIZATION
%==========================================================================
[q_data, epsQ2D_data] = Read_dielectric(eps_file);
[f, Nv_coarse, Nv_dense, N_co, N_de, k] = Read_kmesh(kmesh_file);
N_kpt = size(k, 1); 
[b1, b2] = Reciprocal_lattice_vector(a1, a2);

% --- High-Performance Computing (HPC) Setup ---
poolobj = gcp('nocreate');
if isempty(poolobj) || poolobj.NumWorkers ~= Number_par
    if ~isempty(poolobj), delete(poolobj); end
    parpool('local', Number_par);
end

% Ingest Quantum Coefficients
[Ev, Cv] = READ_VB(bandSource, FunctionPath, Nv, N_kpt, e, NB);
[Ec_kex, Cc_kex] = READ_CB(bandSource, FunctionPath, Nc, N_kpt, e, NB);

%==========================================================================
% [5] BSE HAMILTONIAN (HBSE) CONSTRUCTION
%==========================================================================
[Wq0, wi] = Singularity_and_kweight(N_de, N_co, Nv_coarse, Nv_dense, d, q_data, epsQ2D_data, a, f);
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
        
        parfor j = 1:N_kpt
            % Kinetic Energy Diagonal contribution
            KE_vec = zeros(N_kpt, 1);
            val_KE = (Ec_kex(j, c) - Ev(j, v)) * (eq(c, cc) * eq(v, vv));
            KE_vec(j) = val_KE;
            
            % Momentum-dependent Coulomb Kernel
            q = bsxfun(@minus, k, k(j, :));
            [Vq, epsilon_q] = Screened_Coulomb_interaction(q, b1, b2, d, q_data, epsQ2D_data);
            
            % Direct Interaction integration
            Vd = Direct_Coulomb_interaction(Vq, epsilon_q, k, b1, b2, Wq0(j), wi(j), q, tau, Basis_tau_index, WF_center, N_atoms, Cc_kex, Cv, c, cc, v, vv, j, e, eps0);
            
            HBSE_element(:, j) = KE_vec + (wi(j) / (2 * pi)^2) * Vd;
        end
        
        % Map local blocks to global HBSE
        idx1 = (1:N_kpt) + (((c - 1) + (v - 1) * Nc) * N_kpt);
        idx2 = (1:N_kpt) + (((cc - 1) + (vv - 1) * Nc) * N_kpt);
        HBSE(idx1, idx2) = HBSE_element;
        
        % Symmetrize off-diagonal blocks
       if v ~= vv || c ~= cc
            parfor (j = 1:N_kpt, Number_par)
                HBSE_element(:, j) = HBSE_element(:, j) / (wi(j) / (2 * pi)^2);
            end
            HBSE(idx2, idx1) = (transpose(wi) / (2 * pi)^2) .* ctranspose(HBSE_element);
        end
    end
end
Time_HBSE = toc;
fprintf('[SOLVE] Computation time of constructing BSE matrix: %.2f s\n', Time_HBSE);

%==========================================================================
% [6] EIGENSOLVER (Restored to Original Version)
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
% [7] DATA EXPORT & SYSTEM CLEANUP
%==========================================================================
Total_Time = Time_HBSE + Time_eigs; 
fprintf('[DONE ] Total computation time: %.2f s\n', Total_Time);

save(fullfile(OutputPath, 'A.mat'), 'A', '-v7.3');
save(fullfile(OutputPath, 'Ex.mat'), 'Ex', '-v7.3');
save(fullfile(OutputPath, 'wi.mat'), 'wi', '-v7.3');
k_save = k; k_save(:, 3) = 0; % Project to 2D plane for storage
save(fullfile(OutputPath, 'k.mat'), 'k_save', '-v7.3');

% Generate Performance Logs
output_log(OutputPath, N_kpt, Nv, Nc, Number_par, [Time_HBSE; Time_eigs; Time_HBSE+Time_eigs]);

% Archive and Cleanup 
if isExternalKmesh
    % Archive synchronized bands for research reproducibility
    if userBandsExist
        ArchiveBands = fullfile(OutputPath, 'TB_Bands');
        if ~exist(ArchiveBands, 'dir'), mkdir(ArchiveBands); end
        
        % 1. Move files from Cache to Archive
        movefile(fullfile(CachePath, 'v*_TB_WBSE.txt'), ArchiveBands);
        movefile(fullfile(CachePath, 'c*_TB_WBSE.txt'), ArchiveBands);
        
        % 2. Strip the "_WBSE" suffix to restore original filenames
        archived_files = dir(fullfile(ArchiveBands, '*_WBSE.txt'));
        for f = 1:length(archived_files)
            old_name = archived_files(f).name;
            new_name = strrep(old_name, '_WBSE', '');
            movefile(fullfile(ArchiveBands, old_name), fullfile(ArchiveBands, new_name));
        end
    end
    
    % Final disposal of temporary session files
    temp_cache = fullfile(CachePath, 'kmesh_WBSE.txt');
    if isfile(temp_cache), delete(temp_cache); end
end


% Shutdown parallel environment
delete(gcp('nocreate')); 
fprintf('[DONE ] Execution Successful. Output saved to Exciton_data/.\n');

%==========================================================================
% [AUX] HELPER FUNCTIONS
%==========================================================================
function [fP, aT] = resolve_data_path(fN)
    % Resolves file location prioritizing User-input over Precomputed-cache
    uF = fullfile('User_input', fN);
    if exist(uF, 'file'), fP = uF; aT = 'LOAD '; return; end
    [~, n, e] = fileparts(fN); cF = fullfile('Precomputed_data', [n, '_WBSE', e]);
    if exist(cF, 'file'), fP = cF; aT = 'CACHE'; return; end
    fP = ''; aT = 'GEN  ';
end