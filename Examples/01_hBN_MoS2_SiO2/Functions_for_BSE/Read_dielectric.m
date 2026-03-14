function [q_data, epsQ2D_data] = Read_dielectric(file_path)
% READ_DIELECTRIC Reads q and epsilon data from specific file path.
%
% Input:
%   file_path - Path to the epsilon.txt file
%
% Output:
%   q_data      - Linear momentum column
%   epsQ2D_data - Dielectric function values

    % Handle missing input argument
    if nargin < 1 || isempty(file_path)
        % Fallback to the new standard location
        file_path = 'User_input/epsilon.txt'; 
    end

    % Check if file exists
    if ~isfile(file_path)
        % Second fallback: check the local folder
        if isfile('epsilon.txt')
            file_path = 'epsilon.txt';
        else
            error('Read_dielectric: File not found at %s', file_path);
        end
    end
    
    data = importdata(file_path);
    
    % Extract the first and second columns
    q_data = data(:, 1);
    epsQ2D_data = data(:, 2);
end