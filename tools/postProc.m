% Post processing routine for METiS output 
% POR ENQUANTO TA UMA ZONA. Vou arrumar direitinho depois que validar as
% coisas mais basicas

clear all
% close all
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

flNm = ['..' flSep 'test' flSep 'output_12' flSep 'output.txt'];

flNm_old = '/home/lucas/Google Drive/Doutorado/Software Morison/MATLAB version/OC4_waveOnly_H1000/45deg/firstPlatform_45deg_T10p00_out.mat';
% flNm_old = 'C:\Users\lucas.henrique\Google Drive\Doutorado\Software Morison\MATLAB version\OC4_waveOnly_H1000\45deg\firstPlatform_45deg_T10p00_out.mat';

data = importdata(flNm);
data = data.data;

data_old = load(flNm_old);
time_old = 0:data_old.nump.timeStep:data_old.nump.timeTotal;

figure
set(gcf,'color','w')
sizeOfFont = 12;

subplot(2,3,1)
plot(data(:,1),data(:,2))
hold on
plot(time_old, data_old.strucp_t.pos(:,1))
title('surge')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,2)
plot(data(:,1),data(:,3))
hold on
plot(time_old, data_old.strucp_t.pos(:,2))
title('sway')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,3)
plot(data(:,1),data(:,4))
hold on
plot(time_old, data_old.strucp_t.pos(:,3) - data_old.strucp_t.pos(1,3))
title('heave')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,4)
plot(data(:,1),data(:,5))
hold on
plot(time_old, data_old.strucp_t.pos(:,4))
title('roll')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,5)
plot(data(:,1),data(:,6))
hold on
plot(time_old, data_old.strucp_t.pos(:,5))
title('pitch')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,6)
plot(data(:,1),data(:,7))
hold on
plot(time_old, data_old.strucp_t.pos(:,6))
title('yaw')
set(gca, 'fontsize', sizeOfFont)