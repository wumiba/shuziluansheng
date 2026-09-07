function [V, F] = read_stl_generic(filename)
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open file: %s', filename); end
    header = fread(fid,80,'uint8=>char')';
    frewind(fid);
    is_ascii = startsWith(strtrim(lower(header)), 'solid');
    if is_ascii
        fclose(fid);
        str = fileread(filename);
        vertexTokens = regexp(str, 'vertex\s+([-+eE0-9\.-]+)\s+([-+eE0-9\.-]+)\s+([-+eE0-9\.-]+)', 'tokens');
        verts = cell2mat(cellfun(@(c) str2double(c), vertexTokens, 'UniformOutput', false));
        numVerts = size(verts,1);
        if mod(numVerts,3) ~= 0, error('顶点数不是3的倍数'); end
        rawV = verts;
        rawF = reshape(1:numVerts, 3, [])';
        V = unique(rawV,'rows','stable'); % try preserve unique vertices (simple)
        % simple face rebuilding
        % NOTE: for ASCII STLs returned faces may be rawF. For safety:
        F = rawF;
    else
        fseek(fid,80,'bof');
        numFaces = fread(fid,1,'uint32');
        rawV = zeros(numFaces*3,3);
        rawF = reshape(1:(numFaces*3),3,[])';
        for k=1:numFaces
            fread(fid,3,'float32'); % normal
            v1 = fread(fid,3,'float32')';
            v2 = fread(fid,3,'float32')';
            v3 = fread(fid,3,'float32')';
            rawV((k-1)*3+1,:) = v1;
            rawV((k-1)*3+2,:) = v2;
            rawV((k-1)*3+3,:) = v3;
            fread(fid,1,'uint16'); % attr
        end
        fclose(fid);
        tol = 1e-9;
        [uniqueV,~,idxMap] = uniquetol(rawV, tol, 'ByRows', true);
        V = uniqueV; F = reshape(idxMap, 3, [])';
    end
end