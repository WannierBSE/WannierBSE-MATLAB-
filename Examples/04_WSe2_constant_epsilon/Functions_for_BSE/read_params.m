function params = read_params(fname)
% READ_PARAMS Reads key = value pairs from a text file into a struct.
% Supports:
%   - Scalar numeric values:      Nv = 4
%   - Vectors:                    k_list = 1 2 3
%   - Strings:                    model = monolayer
%   - Inline comments using %     Alpha = 1.5  % fitting param

    fid = fopen(fname, 'r');
    if fid < 0
        error('Cannot open file: %s', fname);
    end

    params = struct();

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        % Skip empty or comment-only lines
        if isempty(line) || startsWith(line, '%')
            continue;
        end

        % Remove inline comment
        commentIdx = strfind(line, '%');
        if ~isempty(commentIdx)
            line = strtrim(line(1:commentIdx(1)-1));
        end

        % Match key = value
        tokens = regexp(line, '^\s*([A-Za-z_]\w*)\s*=\s*(.+)$', 'tokens');
        if isempty(tokens)
            continue;
        end

        key = tokens{1}{1};
        val_str = strtrim(tokens{1}{2});

        % Attempt to parse numeric value or vector
        num_vals = str2num(val_str); %#ok<ST2NM>
        if ~isempty(num_vals)
            params.(key) = num_vals;
        else
            % Store as string if not numeric
            params.(key) = val_str;
        end
    end

    fclose(fid);