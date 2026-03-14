function TB_for_BSE()
    % =========================================================================
    % 1. PARAMETER LOADING (Updated for 1.0.0 Folder Names)
    % =========================================================================
    param_file = '../Parameters/WTB_control.txt';
    params = read_params(param_file);
    VBM = params.VBM; CBM = params.CBM;
    
    main_ctrl = '../Parameters/control.txt';
    main_params = read_params(main_ctrl);
    Nc = main_params.Nc; Nv = main_params.Nv;

    % =========================================================================
    % 2. MODE SELECTION (Updated for 1.0.0 Filenames)
    % =========================================================================
    % We are already inside Precomputed_data/ via WBSE.m
    hr_path = 'wannier90_hr.dat'; 
    
    % Check if ALL required band files exist in the current folder
    files_exist = true;
    for n = 1:Nv
        if ~isfile(['v' num2str(n) '_TB_WBSE.txt']), files_exist = false; break; end
    end
    if files_exist
        for n = 1:Nc
            if ~isfile(['c' num2str(n) '_TB_WBSE.txt']), files_exist = false; break; end
        end
    end
    
    if files_exist
        disp('--- Cache Found: Skipping Band Generation ---');
    elseif isfile(hr_path)
        disp('--- Generating Bands from Hamiltonian ---');
        generate_bands(hr_path, VBM, CBM, Nc, Nv);
    else
        % Check if the file is still in User_input (via the TempInput copy)
        hr_path_alt = '../input_files/wannier90_hr.dat';
        if isfile(hr_path_alt)
             generate_bands(hr_path_alt, VBM, CBM, Nc, Nv);
        else
             error('WannierBSE:NoData', 'wannier90_hr.dat not found.');
        end
    end
end

function generate_bands(hr_path, VBM, CBM, Nc, Nv)
    % 1. Read Structure & K-Mesh (EXACT ORIGINAL LOGIC)
    [a1, a2, a3, ~] = Read_win(); 
    
    kmesh_path = 'kmesh_WBSE.txt'; % Updated filename
    if ~isfile(kmesh_path)
        error('kmesh_WBSE.txt not found.');
    end
    
    % RESTORED: Your original reading method
    k_data = readmatrix(kmesh_path, 'NumHeaderLines', 5);
    k = k_data; 
    Nk = size(k, 1);
    
    % 2. Read Wannier Hamiltonian
    fid = fopen(hr_path, 'r');
    fgetl(fid); 
    NB = str2double(fgetl(fid)); 
    NR = str2double(fgetl(fid));
    gR = fscanf(fid, '%f', NR); 
    hr_data = fscanf(fid, '%f', [7, inf]); 
    hr_data = hr_data'; 
    fclose(fid);
    
    % 3. Construct Hamiltonian
    disp('Constructing and diagonalizing Hamiltonian...');
    R = zeros(NB, NB, NR, 3);
    tij = zeros(NB, NB, NR);
    gR_full = zeros(NB, NB, NR);
    
    for idx = 1:size(hr_data, 1)
        i = hr_data(idx, 4); j = hr_data(idx, 5);
        R_idx = ceil(idx / (NB^2));
        n1 = hr_data(idx, 1); n2 = hr_data(idx, 2); n3 = hr_data(idx, 3);
        R(i, j, R_idx, :) = n1*a1 + n2*a2 + n3*a3;
        tij(i, j, R_idx) = hr_data(idx, 6) + 1i * hr_data(idx, 7);
        gR_full(i, j, R_idx) = gR(R_idx);
    end
    
    % 4. Diagonalization
    E_data = zeros(NB, Nk);
    C_data = zeros(NB, NB, Nk);
    
    tic
    parfor ki = 1:Nk
        % RESTORED: Your exact original reshape
        k_vec = reshape(k(ki, :), [1, 1, 1, 3]); 
        phase = exp(1i * sum(k_vec .* R, 4)); 
        Hk = sum((1 ./ gR_full) .* tij .* phase, 3);
        [C, E] = eig(Hk, 'vector');
        [E_sorted, idx_sort] = sort(real(E));
        E_data(:, ki) = E_sorted;
        C_data(:, :, ki) = C(:, idx_sort);
    end
    t_diag = toc;
    disp(['Diagonalization complete in ' num2str(t_diag) ' s']);
    
    % 5. Write Files (Updated for 1.0.0 suffixes)
    fmt_k = '%20.16f %20.16f %20.16f\n';
    fmt_E = '%20.16f\n';
    fmt_C = '%20.16f %20.16f %20.16f\n';
    
    for n = 1:Nv
        band_idx = VBM - (n - 1);
        fname = ['v' num2str(n) '_TB_WBSE.txt'];
        write_band_file(fname, k, E_data(band_idx, :)', squeeze(C_data(:, band_idx, :)), fmt_k, fmt_E, fmt_C);
    end
    for n = 1:Nc
        band_idx = CBM + (n - 1);
        fname = ['c' num2str(n) '_TB_WBSE.txt'];
        write_band_file(fname, k, E_data(band_idx, :)', squeeze(C_data(:, band_idx, :)), fmt_k, fmt_E, fmt_C);
    end
end

function write_band_file(filename, k, E, C_mat, fmt_k, fmt_E, fmt_C)
    fid = fopen(filename, 'w');
    Nk = size(k, 1);
    nb = size(C_mat, 1);
    for i = 1:Nk
        fprintf(fid, fmt_k, k(i, 1), k(i, 2), k(i, 3));
        fprintf(fid, fmt_E, E(i));
        for orb = 1:nb
            val = C_mat(orb, i);
            fprintf(fid, fmt_C, real(val), imag(val), abs(val));
        end
    end
    fclose(fid);
end