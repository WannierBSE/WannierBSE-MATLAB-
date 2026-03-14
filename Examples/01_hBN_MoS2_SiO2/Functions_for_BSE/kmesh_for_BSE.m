function kmesh_for_BSE()
%KMESH_FOR_BSE Orchestration for k-mesh generation with Deterministic Sorting.
%   This function handles parameter loading, structural data retrieval, 
%   and applies the kmesh_sorter before exporting physical coordinates.

    % 1. INITIALIZATION AND PATH VALIDATION
    output_path = '../Precomputed_data/kmesh_WBSE.txt';
    if isfile(output_path)
        disp(['[INFO] Existing k-mesh file detected: ' output_path]);
        return;
    end

    % 2. CONFIGURATION LOADING
    param_file = '../Parameters/kmesh_control.txt';
    if ~isfile(param_file)
        error('WannierBSE:FileNotFound', 'Configuration file not found: %s', param_file);
    end
    
    params = read_params(param_file);
    Nk_coarse   = params.Nk_coarse;
    Nk_dense    = params.Nk_dense;
    Dense_range = params.Dense_range;
    
    if Dense_range > Nk_coarse / 2
        error('WannierBSE:ValidationFailed', ...
              'Dense_range must be <= Nk_coarse/2 to prevent point duplication.');
    end

    % 3. LATTICE DATA ACQUISITION & ENGINE EXECUTION
    [a1, a2, ~, a_const] = Read_win(); 
    
    % Execute core geometric engine to obtain dimensionless coordinates
    [final_kmesh, total_coarse, total_dense, ~] = ...
        kmesh_hexagonal_engine(Nk_coarse, Nk_dense, Dense_range, a1, a2, a_const);

    % =========================================================================
    % NEW: INTEGRATION OF DETERMINISTIC SORTING
    % =========================================================================
    % Points are sorted before scaling to ensure dimensionless precision.
    [final_kmesh, ~] = kmesh_sorter(final_kmesh, total_dense, total_coarse);

    % 4. EXPORT TO KMESH.TXT (Physical Units Application)
    fid_out = fopen(output_path, 'w');
    if fid_out == -1
        error('WannierBSE:FileWriteError', 'Unable to open %s for writing.', output_path);
    end

    % Write Header Information
    fprintf(fid_out, '%d\n', Dense_range);
    fprintf(fid_out, '%d\n', Nk_coarse);
    fprintf(fid_out, '%d\n', Nk_dense);
    fprintf(fid_out, '%d\n', total_coarse); 
    fprintf(fid_out, '%d\n', total_dense);

    % Scale factor for physical K-space conversion (1/Angstrom)
    scale_factor = pi / a_const;

    % Export 3D k-points with high-fidelity precision
    for i = 1:size(final_kmesh, 1)
        kx_scaled = final_kmesh(i,1) * scale_factor;
        ky_scaled = final_kmesh(i,2) * scale_factor;
        fprintf(fid_out, '%20.16f    %20.16f     %20.16f\n', kx_scaled, ky_scaled, 0.0);
    end
    fclose(fid_out);
end