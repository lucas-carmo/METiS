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


% flNm = '/home/lucas/Nextcloud/Doutorado/Casos teste/preliminares/OC4-regular/HD_movel-HS_linear/H2p0m_T30p00_out.txt';
flNm = '..\test\OC4_200m_45deg_Tirreg_out.txt';
flNm2 = '..\test\OC4_200m_45deg_Tirreg_out.txt';



legStr = {'HD = 3', 'HD = 2'};



%=========================================================================%
data_root = importdata(flNm);
data = data_root.data;

data_2_root = importdata(flNm2);
data_2 = data_2_root.data;

figure
set(gcf,'color','w')
sizeOfFont = 12;

subplot(2,3,1)
plot(data(:,1),data(:,3))
hold on
plot(data_2(:,1),data_2(:,3))
hl = legend(legStr);
set(hl, 'location', 'northwest')
title('surge')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,2)
plot(data(:,1),data(:,4))
hold on
plot(data_2(:,1),data_2(:,4))
title('sway')   
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,3)
plot(data(:,1),data(:,5))
hold on
plot(data_2(:,1),data_2(:,5))
title('heave')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,4)
plot(data(:,1),data(:,6))
hold on
plot(data_2(:,1),data_2(:,6))
title('roll')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,5)
plot(data(:,1),data(:,7))
hold on
plot(data_2(:,1),data_2(:,7))
title('pitch')
set(gca, 'fontsize', sizeOfFont)

subplot(2,3,6)
plot(data(:,1),data(:,8))
hold on
plot(data_2(:,1),data_2(:,8))
title('yaw')
set(gca, 'fontsize', sizeOfFont)


figure
plot(data(:,1),data(:,2))
ylabel('wave elevation (m)')
xlabel('time (s)')


% % Hydrodynamic force - Total
% figure
% set(gcf,'color','w')
% sizeOfFont = 12;
% 
% subplot(2,3,1)
% plot(data(:,1),data(:,20))
% hold on
% plot(data_2(:,1),data_2(:,20))
% legend(legStr)
% title('force-hd-surge')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,2)
% plot(data(:,1),data(:,21))
% hold on
% plot(data_2(:,1),data_2(:,21))
% title('force-hd-sway')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,3)
% plot(data(:,1),data(:,22))
% hold on
% plot(data_2(:,1),data_2(:,22))
% title('force-hd-heave')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,4)
% plot(data(:,1),data(:,23))
% hold on
% plot(data_2(:,1),data_2(:,23))
% title('force-hd-roll')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,5)
% plot(data(:,1),data(:,24))
% hold on
% plot(data_2(:,1),data_2(:,24))
% title('force-hd-pitch')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,6)
% plot(data(:,1),data(:,25))
% hold on
% plot(data_2(:,1),data_2(:,25))
% title('force-hd-yaw')
% set(gca, 'fontsize', sizeOfFont)
% 
% 
% 
% % Hydrodynamic force - Inertia
% figure
% set(gcf,'color','w')
% sizeOfFont = 12;
% 
% subplot(2,3,1)
% plot(data(:,1),data(:,26))
% hold on
% plot(data_2(:,1),data_2(:,26), '--')
% legend(legStr)
% title('force-hd-inertia-surge')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,2)
% plot(data(:,1),data(:,27))
% hold on
% plot(data_2(:,1),data_2(:,27), '--')
% title('force-hd-inertia-sway')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,3)
% plot(data(:,1),data(:,28))
% hold on
% plot(data_2(:,1),data_2(:,28), '--')
% title('force-hd-inertia-heave')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,4)
% plot(data(:,1),data(:,29))
% hold on
% plot(data_2(:,1),data_2(:,29), '--')
% title('force-hd-inertia-roll')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,5)
% plot(data(:,1),data(:,30))
% hold on
% plot(data_2(:,1),data_2(:,30), '--')
% title('force-hd-inertia-pitch')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,6)
% plot(data(:,1),data(:,31))
% hold on
% plot(data_2(:,1),data_2(:,31), '--')
% title('force-hd-inertia-yaw')
% set(gca, 'fontsize', sizeOfFont)
% 
% 
% % Hydrodynamic force - Drag
% figure
% set(gcf,'color','w')
% sizeOfFont = 12;
% 
% subplot(2,3,1)
% plot(data(:,1),data(:,32))
% hold on
% plot(data_2(:,1),data_2(:,32), '--')
% legend(legStr)
% title('force-hd-drag-surge')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,2)
% plot(data(:,1),data(:,33))
% hold on
% plot(data_2(:,1),data_2(:,33), '--')
% title('force-hd-drag-sway')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,3)
% plot(data(:,1),data(:,34))
% hold on
% plot(data_2(:,1),data_2(:,34), '--')
% title('force-hd-drag-heave')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,4)
% plot(data(:,1),data(:,35))
% hold on
% plot(data_2(:,1),data_2(:,35), '--')
% title('force-hd-drag-roll')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,5)
% plot(data(:,1),data(:,36))
% hold on
% plot(data_2(:,1),data_2(:,36), '--')
% title('force-hd-drag-pitch')
% set(gca, 'fontsize', sizeOfFont)
% 
% subplot(2,3,6)
% plot(data(:,1),data(:,37))
% hold on
% plot(data_2(:,1),data_2(:,37), '--')
% title('force-hd-drag-yaw')
% set(gca, 'fontsize', sizeOfFont)