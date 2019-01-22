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

flNm = '/home/lucas/Nextcloud/Tese/Casos teste/OC4/metis_OC4_completo_h200_45deg/OC4_200m_45deg_T30p00_out.txt';
flNm2 = '/home/lucas/Nextcloud/Tese/Casos teste/OC4/metis_OC4_completo_h200_45deg_fixed/OC4_200m_45deg_fixed_T30p00_out.txt';


legStr = {'Livre', 'Fixo'};

data = importdata(flNm);
data = data.data;

data_2 = importdata(flNm2);
data_2 = data_2.data;


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


% figure
% plot(data(:,1), data(:,24))
% hold on
% plot(time_old, data_old.moment(:,1))
% title('force total roll')


% figure
% GM = 10.67;
% g = 9.81;
% delta = 1.43e7 * g;
% 
% plot(data(:,1), data(:,11), 'b')
% hold on
% plot(time_old, data_old.momentHydrostatic(:,1), 'r')
% hold on
% plot(data(1:20:end,1), -delta * GM * sin(data(1:20:end,5)), '--ob')
% hold on
% plot(time_old(1:20:end), -delta * GM * sin(data_old.strucp_t.pos(1:20:end,4)), '--or')
% 
% title('hydrostatic roll')
% 
% 
% figure
% % Awl = pi*6.5^2/4 + 3*pi*12^2/4;
% % rho = 1025;
% K_heave = 3.84e6; % N/m
% 
% plot(data(:,1), data(:,10) - data(1,10), 'b')
% hold on
% plot(time_old(1:end-1), data_old.forceHydrostatic(1:end-1,3) - data_old.forceHydrostatic(1,3), 'r')
% hold on
% plot(data(1:20:end,1), - K_heave * data(1:20:end,4), '--ob')
% hold on
% plot(time_old(1:20:end-1), - K_heave * (data_old.strucp_t.pos(1:20:end-1,3) - data_old.strucp_t.pos(1,3)), '--or')
% title('hydrostatic heave')


% figure
% plot(data(:,1), data(:,23))
% hold on
% plot(time_old, data_old.momentHeavePlate(:,1))
% title('moment heave plate roll')
% 
% figure
% plot(data(:,1), data(:,22))
% hold on
% plot(time_old, data_old.forceHeavePlate(:,3))
% title('force heave plate heave')