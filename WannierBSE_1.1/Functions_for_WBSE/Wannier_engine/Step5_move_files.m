%==========================================================================
% STEP5_MOVE_FILES
%
% DESCRIPTION:
%   Moves normalized Wannier-function data from the temporary workspace to the permanent WBSE precomputed-data cache.
%==========================================================================

load(fullfile(Paths.temp_normalized, 'W_up.mat'))
load(fullfile(Paths.temp_normalized, 'W_down.mat'))
load(fullfile(Paths.temp_normalized, 'r.mat'))

% Create final destination folder if it doesn't exist
if ~exist(Paths.final_output, 'dir')
    mkdir(Paths.final_output);
end

% Save results to the final permanent directory (Precomputed_data/Wannier_functions_WBSE/)
save(fullfile(Paths.final_output, 'WF_up.mat'),'W_up','-v7.3')
save(fullfile(Paths.final_output, 'WF_down.mat'),'W_down','-v7.3')
save(fullfile(Paths.final_output, 'r.mat'),'r','-v7.3')

% Copy metadata files to final destination
copyfile(fullfile(Paths.temp_normalized, 'GridInfo.txt'), fullfile(Paths.final_output, 'GridInfo.txt'))
copyfile(fullfile(Paths.temp_normalized, 'WF_AtomPosition.txt'), fullfile(Paths.final_output, 'WF_AtomPosition.txt'))

%-------delete data in previous step (Normalized_WFs)--------
if exist(fullfile(Paths.temp_normalized, 'W_up.mat'), 'file'), delete(fullfile(Paths.temp_normalized, 'W_up.mat')); end
if exist(fullfile(Paths.temp_normalized, 'W_down.mat'), 'file'), delete(fullfile(Paths.temp_normalized, 'W_down.mat')); end
if exist(fullfile(Paths.temp_normalized, 'r.mat'), 'file'), delete(fullfile(Paths.temp_normalized, 'r.mat')); end

% Final Workspace Cleanup: Remove the entire temporary directory
% This realizes the "zero-modification" requirement by wiping the scratch area.
if exist(Paths.workspace, 'dir')
    rmdir(Paths.workspace, 's');
    fprintf('Temporary workspace cleaned up successfully.\n');
end
