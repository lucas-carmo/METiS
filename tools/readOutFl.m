% Read METiS results (_out file)
%
% out = readOutFl(flNm) reads a METiS _out file with the results of a time
% domain simulation. The file must consist of columns with headers in the
% first line and the data in the rest of the file:
%
%
%        Time       Surge              Sway             Heave              Roll             Pitch               Yaw
%   0.0000e+00    0.0000e+00        0.0000e+00        0.0000e+00        0.0000e+00        0.0000e+00        0.0000e+00
%   2.5000e-01    1.0000e-08        1.0000e-08        1.0000e-08        1.0000e-08        1.0000e-08        1.0000e-08

function out = readOutFl(flNm)

if ~ischar(flNm)
    error('The input must be a string with the path to a single METiS out file')
end


data = importdata(flNm);

if ~isfield(data, 'colheaders')
    error('Invalid file format. Data must be arranged in columns with headers in the first line and the data in the rest of the file')
end


for ii = 1 : numel(data.colheaders)
    str4field = data.colheaders{ii};
    out.(str4field) = data.data(:,ii);
end