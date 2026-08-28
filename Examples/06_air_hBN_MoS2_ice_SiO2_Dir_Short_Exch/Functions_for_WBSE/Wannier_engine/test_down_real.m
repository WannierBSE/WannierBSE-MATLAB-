%==========================================================================
% TEST_DOWN_REAL
%
% DESCRIPTION:
%   Regression helper for parsing and exporting the spin-down real XSF component.
%==========================================================================

xsf_downreal_path = fullfile(Paths.xsf_root, 'down_real');
output_target_dir = fullfile(Paths.temp_output, 'down_real');

% --- Verify Input Data ---
% Validate the existence of the source component folder
if ~isfolder(xsf_downreal_path)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', xsf_downreal_path);
    error(errorMessage);
end

% --- Directory Setup ---
% Ensure the specific output container exists within the workspace
if ~exist(output_target_dir, 'dir')
    mkdir(output_target_dir);
end

% --- File Discovery & Sequencing ---
% Retrieve and sort all .xsf files to ensure consistent Wannier function indexing
filePattern1 = fullfile(xsf_downreal_path, '*.xsf');
theFiles1 = dir(filePattern1);
baseFileNames = natsortfiles({theFiles1.name});
LEN = length(theFiles1);

% --- Primary Extraction Loop ---
for k = 1:LEN
    filename = baseFileNames{k};
    fullFilePath = fullfile(xsf_downreal_path, filename);
    fileID = fopen(fullFilePath, 'r');
    
    % Locate data grid boundaries within the XSF format
    A = regexp(fileread(fullFilePath), '\n', 'split');
    end_grid_row = find(contains(A, 'END_DATAGRID_3D'));
    data_dim_row = find(contains(A, 'DATAGRID_3D_UNKNOWN'));
    delimiter = ' ';
    
    % Define data extraction range
    startRow = data_dim_row + 6;
    endRow = end_grid_row - 1;
    formatSpec = '%f%f%f%f%f%f%[^\n\r]';
    
    % Read volumetric data
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, ...
        'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
        'TextType', 'string', 'HeaderLines', startRow-1, ...
        'ReturnOnError', false, 'EndOfLine', '\r\n');
    raw = [dataArray{1:end-1}];
    fclose(fileID);
    
    % Extract grid dimensions (nx, ny, nz)
    fileID = fopen(fullFilePath, 'r');
    data_dim = textscan(fileID, '%f %f %f', 1, 'Delimiter', delimiter, ...
        'MultipleDelimsAsOne', true, 'headerlines', data_dim_row);
    fclose(fileID);
    
    % Extract unit cell geometry
    fileID = fopen(fullFilePath, 'r');
    cell0 = textscan(fileID, '%f %f %f', 1, 'headerlines', data_dim_row+1);
    fclose(fileID);
    fileID = fopen(fullFilePath, 'r');
    cell1 = textscan(fileID, '%f %f %f', 1, 'headerlines', data_dim_row+2);
    fclose(fileID);
    fileID = fopen(fullFilePath, 'r');
    cell2 = textscan(fileID, '%f %f %f', 1, 'headerlines', data_dim_row+3);
    fclose(fileID);
    fileID = fopen(fullFilePath, 'r');
    cell3 = textscan(fileID, '%f %f %f', 1, 'headerlines', data_dim_row+4);
    fclose(fileID);
    
    A0 = cell2mat(cell0);
    A1 = cell2mat(cell1);
    A2 = cell2mat(cell2);
    A3 = cell2mat(cell3);
    
    % --- Data Reshaping & Normalization ---
    nx = data_dim{1,1};
    ny = data_dim{1,2};
    nz = data_dim{1,3};
    
    temp = raw';
    data = temp(:);
    
    % Create the coordinate and value matrix (position1)
    position1 = zeros(nx*ny*nz, 7);
    x = (1:nx)';
    x_rep = repmat(x, ny*nz, 1);
    
    y = (1:ny)';
    y_rep = repelem(y, nx);
    y_rep = repmat(y_rep, nz, 1);
    
    z = (1:nz)';
    z_rep = repelem(z, nx*ny);
    
    for n = 1:nx*ny*nz
        position1(n,1) = x_rep(n);
        position1(n,2) = y_rep(n);
        position1(n,3) = z_rep(n);
        
        % Map grid indices to Cartesian coordinates
        position1(n,4) = A1(1)*(x_rep(n)-1)/(nx-1) + A2(1)*(y_rep(n)-1)/(ny-1) + A3(1)*(z_rep(n)-1)/(nz-1) + A0(1);
        position1(n,5) = A1(2)*(x_rep(n)-1)/(nx-1) + A2(2)*(y_rep(n)-1)/(ny-1) + A3(2)*(z_rep(n)-1)/(nz-1) + A0(2);
        position1(n,6) = A1(3)*(x_rep(n)-1)/(nx-1) + A2(3)*(y_rep(n)-1)/(ny-1) + A3(3)*(z_rep(n)-1)/(nz-1) + A0(3);
        position1(n,7) = data(n,1);
    end
    
    % --- Save Result ---
    % Export processed matrix to the temporary workspace
    save(fullfile(output_target_dir, [num2str(k), '_WF.mat']), 'position1', '-v7.3');
    
    %--------save dat file
%     fid = fopen([num2str(k),'_WF.dat'],'w');
%     for i=1:size(position1,1)
%         fprintf(fid,'%4i %4i %4i %11.6f %11.6f %11.6f %11.10f',...
%            position1(i,1),position1(i,2),position1(i,3),position1(i,4),position1(i,5),position1(i,6),position1(i,7));
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
    %--------save dat file
    
    clearvars raw data position1 x_rep y_rep z_rep temp
end
