function raw_write(V, filename)
% RAW_WRITE(V, filename), write the image V to a RAW file (binary).
% The datatype will be inferred from V, and the endianness will be
% your machine's default.
%
% RAW_WRITE will append size, datatype and endianness to the filename,
% for example: 'filename_100x200x300_uint8_l.raw'.
% This information is often required when using a RAW file, it is also
% interpreted by raw_read.

% get machine format
[~, ~, endianness] = computer;
endianness = lower(endianness);

% get the datatype
datatype = class(V);

% split filename
[pathname, basename, ext] = fileparts(filename);

% add RAW extension
if isempty(ext)
    ext = '.raw';
end

% append meta-data
basename = sprintf('%s_%s_%s_%dx%dx%d', basename, endianness, datatype, size(V));

% update file name
filename = fullfile(pathname, sprintf('%s%s', basename, ext));

% write to file
fid = fopen(filename, 'w+');
fwrite(fid, V, datatype, 0, endianness);
fclose(fid);

end
