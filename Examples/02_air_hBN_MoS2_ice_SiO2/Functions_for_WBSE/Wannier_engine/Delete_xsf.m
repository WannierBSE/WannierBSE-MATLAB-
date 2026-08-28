%==========================================================================
% DELETE_XSF
%
% DESCRIPTION:
%   Deletes intermediate Wannier XSF component files after conversion while preserving the first reference XSF in each spin component folder.
%==========================================================================

if exist('cleanup_xsf', 'var') && cleanup_xsf == false
    fprintf('  -> XSF cleanup skipped; source .xsf files retained.\n');
else
    fprintf('  -> Removing redundant intermediate .xsf files...\n');
    
    % --- Safety Check ---
    if ~exist('Paths', 'var')
        error('Delete_xsf: Path structure not found. Please run init_wannier_paths first.');
    end

    % Define component subdirectories
    subfolders = {'up_real', 'up_imag', 'down_real', 'down_imag'};
    
    for xi = 1:numel(subfolders)
        target_dir = fullfile(Paths.xsf_root, subfolders{xi});
        
        if isfolder(target_dir)
            % Identify all relevant Wannier files in the current component folder
            all_files = dir(fullfile(target_dir, 'wannier90_*.xsf'));
            
            % --- Targeted Deletion ---
            % We iterate from the 2nd file onward to preserve 'wannier90_00001.xsf'
            % as a reference for visualization or metadata verification.
            for i = 2:length(all_files)
                % Generate the 5-digit zero-padded filename (e.g., wannier90_00002.xsf)
                filename = sprintf('wannier90_%05d.xsf', i);
                filepath = fullfile(target_dir, filename);
                
                if exist(filepath, 'file')
                    delete(filepath);
                end
            end
        end
    end
    fprintf('  -> Intermediate .xsf files removed from User_input.\n');
end
