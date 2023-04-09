function V = raw_read(filename, im_size, datatype, endianness)
% V = RAW_READ(filename, im_size, datatype, endianness), reads RAW file.
% The input 'im_size' is a 3-element vector that defines the matrix size.
% See raw_write for details on 'datatype' and 'endianness' parameters.
%
% V = RAW_READ(filename, im_size, datatype), uses 'native' endianness.
%
% V = RAW_READ(filename), parses the filename, uses the same format as
% raw_write, ex: 'filename_200x200x200_uint8_n.raw'

if nargin < 2

    % split filename
    [~, basename, ~] = fileparts(filename);

    % parse filename
    a = regexp(basename, '(.*)_(.*)_(.*)_(\d*)x(\d*)x(\d*)', 'tokens');
    assert(not(isempty(a)), ...
        'Error interpreting the file name, please provide the meta-data');

    a = a{1};
    endianness = a{2};
    datatype = a{3};
    im_size = [str2double(a{4}), str2double(a{5}), str2double(a{6})];

end

if nargin < 4
    endianness = 'n';
end

if any(strcmpi(datatype, {'int8', 'uint8', 'char'}))
    word = 1;
elseif any(strcmpi(datatype, {'int16', 'uint16'}))
    word = 2;
elseif any(strcmpi(datatype, {'int32', 'uint32', 'single', 'float32'}))
    word = 4;
elseif any(strcmpi(datatype, {'int64', 'uint64', 'double', 'float64'}))
    word = 8;
else
    error('unknown data type');
end

fid = fopen(filename, 'r');
fseek(fid, 0, 'eof'); % go to the end of the file
nBytes = ftell(fid);
fseek(fid, 0, 'bof'); % rewind to start of the file

assert(prod(im_size) * word == nBytes, ...
    'File size does not match specified meta-data.');

% read the data
V = fread(fid, inf, ['*' datatype], 0, endianness);
fclose(fid);

% shape it
V = reshape(V, im_size);
