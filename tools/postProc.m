% Post processing routine for METiS output 
% POR ENQUANTO TA UMA ZONA. Vou arrumar direitinho depois que validar as
% coisas mais basicas

clear all
close all
clc

% I do not know why, but MATLAB's function 'filesep' was not working
% properly in Ubuntu, so I decided to to specify the file separator this
% way
if ispc
    flSep = '\';
elseif isunix
    flSep = '/';
elseif ismac
    flSep = ':';
else
    flSep = '/';
end

flNm = ['..' flSep 'test' flSep 'output_1' flSep 'output.txt'];

flNm_old = ['/home/lucas/Google Drive/Doutorado/Software Morison/MATLAB version/OC4_waveOnly_H1000/45deg/firstPlatform_45deg_T10p00_out.mat'];

data = importdata(flNm);
data = data.data;

data_old = load(flNm_old);

figure
plot(data(:,1),data(:,2))




