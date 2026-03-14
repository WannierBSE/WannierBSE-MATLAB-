function [WF_center, Basis_tau_index_cell] = Read_WF_center(N_atoms) % Modified by Oscar
%
% This version is modified to return a cell array
% which is compatible with the Python-MATLAB engine.
%
% find number of rows
fid = fopen('../parameters/WF_centers.txt','r');
Nh = 0;
while fgetl(fid)~=-1
    Nh = Nh+1;
end
frewind(fid)
WF_center_data = strings(Nh,1);
for i=1:Nh
WF_center_data{i} = fgetl(fid);
end
for i=1:(size(WF_center_data,1)-2)
WF_center(i,:) = str2num(WF_center_data{i+2});
end
fclose(fid);

% --- MODIFICATION ---
% Create a cell array, which Python can understand as a list
Basis_tau_index_cell = cell(N_atoms, 1); % Modified by Oscar
for Ni=1:N_atoms
    index_tau_i = find(WF_center(:,2)==Ni);
    Basis_tau_index_cell{Ni} = index_tau_i; % Store in a cell {Ni} % Modified by Oscar
end
% --- END MODIFICATION ---
end