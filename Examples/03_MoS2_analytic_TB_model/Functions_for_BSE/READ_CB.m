function [Ec, Cc] = READ_CB(bandSource, FunctionPath, Nc, N_kpt, e, NB)
    Ec = zeros(N_kpt, Nc);
    Cc = zeros(N_kpt, NB, Nc);
    
    for i = 1:Nc
        CB_NAME = fullfile(bandSource, sprintf('c%d_TB.txt', i));
        if ~isfile(CB_NAME)
            CB_NAME = fullfile(bandSource, sprintf('c%d_TB_WBSE.txt', i));
        end
        
        CB_data_raw = importdata(CB_NAME);
        if isstruct(CB_data_raw), CB_data_raw = CB_data_raw.data; end
        
        [Ec(:,i), Cc(:,:,i)] = Read_TBdata(N_kpt, CB_data_raw, e, NB);
    end
end