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


% flNm = '/home/lucas/Nextcloud/Tese/Casos teste/OC4/metis_OC4_completo_h200_45deg/OC4_200m_45deg_T30p00_out.txt';
% flNm2 = '/home/lucas/Nextcloud/Tese/Casos teste/OC4/metis_OC4_completo_h200_45deg_fixed/OC4_200m_45deg_fixed_T30p00_out.txt';

flNm = 'C:\Users\lucas.henrique\Nextcloud\Tese\Casos teste\OC4\metis_OC4_completo_h200_45deg\OC4_200m_45deg_T30p00_out.txt';
flNm2 = 'C:\Users\lucas.henrique\Nextcloud\Tese\Casos teste\OC4\metis_OC4_completo_h200_45deg_fixed\OC4_200m_45deg_fixed_T30p00_out.txt';

legStr = {'Livre', 'Fixo'};



%=========================================================================%
data_root = importdata(flNm);
data = data_root.data;

data_2_root = importdata(flNm2);
data_2 = data_2_root.data;


set(gcf,'color','w')
sizeOfFont = 12;

subplot(2,3,1)
plot(data(:,1),data(:,2))
hold on
plot(data_2(:,1),data_2(:,2))
legend(legStr)
title('surge')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,2)
plot(data(:,1),data(:,3))
hold on
plot(data_2(:,1),data_2(:,3))
title('sway')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,3)
plot(data(:,1),data(:,4))
hold on
plot(data_2(:,1),data_2(:,4))
title('heave')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,4)
plot(data(:,1),data(:,5))
hold on
plot(data_2(:,1),data_2(:,5))
title('roll')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,5)
plot(data(:,1),data(:,6))
hold on
plot(data_2(:,1),data_2(:,6))
title('pitch')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,6)
plot(data(:,1),data(:,7))
hold on
plot(data_2(:,1),data_2(:,7))
title('yaw')
set(gca, 'fontsize', sizeOfFont)





% Forces
figure
set(gcf,'color','w')
sizeOfFont = 12;

subplot(2,3,1)
plot(data(:,1),data(:,20))
hold on
plot(data_2(:,1),data_2(:,20))
legend(legStr)
title('force-hd-surge')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,2)
plot(data(:,1),data(:,21))
hold on
plot(data_2(:,1),data_2(:,21))
title('force-hd-sway')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,3)
plot(data(:,1),data(:,22))
hold on
plot(data_2(:,1),data_2(:,22))
title('force-hd-heave')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,4)
plot(data(:,1),data(:,23))
hold on
plot(data_2(:,1),data_2(:,23))
title('force-hd-roll')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,5)
plot(data(:,1),data(:,24))
hold on
plot(data_2(:,1),data_2(:,24))
title('force-hd-pitch')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,6)
plot(data(:,1),data(:,25))
hold on
plot(data_2(:,1),data_2(:,25))
title('force-hd-yaw')
set(gca, 'fontsize', sizeOfFont)