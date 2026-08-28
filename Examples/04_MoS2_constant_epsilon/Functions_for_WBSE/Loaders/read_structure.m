function [a1, a2, a3, a, tau, N_atoms] = read_structure(varargin)
%==========================================================================
% READ_STRUCTURE
%
% DESCRIPTION:
%   Loads lattice vectors and atomic coordinates from Wannier90 or legacy structure input files.
%
% INPUT ARGUMENTS:
%   varargin : optional User_input directory path.
%
% OUTPUT:
%   a1, a2, a3 : real-space lattice vectors.
%   a : in-plane lattice constant norm.
%   tau : atomic Cartesian coordinates.
%   N_atoms : number of atoms in the unit cell.
%==========================================================================

persistent structure_source_reported

if nargin > 0
    user_input_path = varargin{1};
else
    loader_path = fileparts(mfilename('fullpath'));
    project_root = fullfile(loader_path, '..', '..');
    user_input_path = fullfile(project_root, 'User_input');
end

wannier90_path = fullfile(user_input_path, 'wannier90.win');
structure_path = fullfile(user_input_path, 'structure.txt');

if isfile(wannier90_path)
    source_path = wannier90_path;
    source_name = fullfile('User_input', 'wannier90.win');
elseif isfile(structure_path)
    source_path = structure_path;
    source_name = fullfile('User_input', 'structure.txt');
else
    error('WannierBSE:StructureNotFound', ...
        'No structure input was found. Expected either %s or %s.', ...
        wannier90_path, structure_path);
end

fid = fopen(source_path, 'r');
if fid == -1
    error('WannierBSE:StructureOpenFailed', ...
        'Cannot open structure input: %s', source_path);
end
cleanup_obj = onCleanup(@() fclose(fid));

lines = strings(0, 1);
while ~feof(fid)
    lines(end + 1, 1) = string(fgetl(fid)); %#ok<AGROW>
end

lattice_begin = find(strtrim(lines) == "begin unit_cell_cart", 1);
atoms_begin = find(strtrim(lines) == "begin atoms_cart", 1);
atoms_end = find(strtrim(lines) == "end atoms_cart", 1);

if isempty(lattice_begin) || lattice_begin + 3 > numel(lines)
    error('WannierBSE:InvalidStructure', ...
        'Missing or incomplete unit_cell_cart block in %s.', source_path);
end
if isempty(atoms_begin) || isempty(atoms_end) || atoms_end <= atoms_begin + 1
    error('WannierBSE:InvalidStructure', ...
        'Missing or empty atoms_cart block in %s.', source_path);
end

a1 = parse_vector(lines(lattice_begin + 1), source_path);
a2 = parse_vector(lines(lattice_begin + 2), source_path);
a3 = parse_vector(lines(lattice_begin + 3), source_path);
a = norm(a1);

N_atoms = atoms_end - atoms_begin - 1;
tau = zeros(N_atoms, 3);
for atom_index = 1:N_atoms
    atom_values = sscanf(char(lines(atoms_begin + atom_index)), '%*s %f %f %f');
    if numel(atom_values) ~= 3
        error('WannierBSE:InvalidStructure', ...
            'Invalid atomic coordinate in %s at atom %d.', ...
            source_path, atom_index);
    end
    tau(atom_index, :) = atom_values.';
end

if isempty(structure_source_reported)
    fprintf('[INFO ] Structure extracted from %s.\n', source_name);
    structure_source_reported = true;
end

clear cleanup_obj
end


function vector = parse_vector(line, source_path)
%==========================================================================
% PARSE_VECTOR
%
% DESCRIPTION:
%   Parses one lattice-vector line from a structure file.
%
% OUTPUT:
%   vector : three-component row vector.
%==========================================================================

vector = sscanf(char(line), '%f %f %f');
if numel(vector) ~= 3
    error('WannierBSE:InvalidStructure', ...
        'Invalid lattice vector in %s.', source_path);
end
vector = vector.';
end
