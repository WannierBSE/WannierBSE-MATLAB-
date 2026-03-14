function [a1,a2,a3,a,tau,N_atoms] = Read_win
% Open file
fid = fopen('../Parameters/structure.txt','r');
if fid == -1
    error('Cannot open file "structure".');
end

% Count the number of lines
Nh = 0;
while ~feof(fid)
    fgetl(fid);
    Nh = Nh + 1;
end
frewind(fid); % Go back to start of file


% Read all lines into cell array
win_data = strings(Nh,1);
for i=1:Nh
    win_data{i} = fgetl(fid);
end
fclose(fid);


Line_lattice_vectors_begin = find(win_data=='begin unit_cell_cart');
a1 = str2num(win_data{Line_lattice_vectors_begin+1});
a2 = str2num(win_data{Line_lattice_vectors_begin+2});
a3 = str2num(win_data{Line_lattice_vectors_begin+3});

a = norm(a1);

Line_atomes_begin = find(win_data=='begin atoms_cart');
Line_atomes_end = find(win_data=='end atoms_cart');
N_atoms = Line_atomes_end-Line_atomes_begin-1;

tau = zeros(N_atoms,3);
for i=1:N_atoms
tau(i,:) = str2num(win_data{Line_atomes_begin+i}(4:end));
end

end