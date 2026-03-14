function [Ev, Cv] = READ_VB(bandSource, FunctionPath, Nv, N_kpt, e, NB)
    % Pre-allocate to prevent size errors
    Ev = zeros(N_kpt, Nv);
    Cv = zeros(N_kpt, NB, Nv);
    
    for i = 1:Nv
        % Use the 1.0.0 Naming Convention
        VB_NAME = fullfile(bandSource, sprintf('v%d_TB.txt', i));
        
        % Fallback for Cache naming
        if ~isfile(VB_NAME)
            VB_NAME = fullfile(bandSource, sprintf('v%d_TB_WBSE.txt', i));
        end
        
        % Load the raw numeric data
        VB_data_raw = importdata(VB_NAME);
        if isstruct(VB_data_raw), VB_data_raw = VB_data_raw.data; end
        
        % USE YOUR ORIGINAL HELPER FUNCTION
        % This function handles the complex number conversion and column splitting
        [Ev(:,i), Cv(:,:,i)] = Read_TBdata(N_kpt, VB_data_raw, e, NB);
    end
end