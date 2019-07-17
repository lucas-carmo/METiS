% Post processing routine for METiS output 

clear all
close all
clc


flNm = {'C:\Users\lucas.henrique\Nextcloud\Doutorado\Casos teste\preliminares\OC4-regular\HD_movel-HS_linear\H2p0m_T30p00_out.txt'
        'C:\Users\lucas.henrique\Nextcloud\Doutorado\Casos teste\preliminares\OC4-regular\HD_movel-HS_linear\H4p0m_T30p00_out.txt'
        'C:\Users\lucas.henrique\Nextcloud\Doutorado\Casos teste\preliminares\OC4-regular\HD_movel-HS_linear\H8p0m_T30p00_out.txt'
        };

legSpec = { 'H = 2m'
            'H = 4m'
            'H = 8m'
          };
    
%===== Choose the output
activeDoFs = [1 1 1 1 1 1];
fowt_disp = 1;
fowt_vel = 1;
fowt_acc = 1;

hd_inertia_force = 1;
hd_drag_force = 1;
hd_fk_force = 1;
hd_force = 1;
hs_force = 1;
ad_hub_force = 1;
total_force = 1;

% Wave outputs
wave_elev = 1;
wave_vel = 1;
wave_acc = 1;
wave_press = 1;



%===== Plot style
width4Line = 2;
colors4Plot = num2cell(get(groot,'defaultAxesColorOrder'), 2);
sizeOfFont = 12;


%=========================================================================%
totalDoF = sum(activeDoFs);

for ii = 1:numel(flNm)
    data = readOutFl(flNm{ii});
    
    
    if fowt_disp == 1
        if ii == 1
            fig_fowt_disp = figure('name', 'FOWT displacement' ,'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);                       
            ax_fowt_disp_surge = subplot(2,3,1);    
            ax_fowt_disp_sway = subplot(2,3,2);    
            ax_fowt_disp_heave = subplot(2,3,3);    
            ax_fowt_disp_roll = subplot(2,3,4);    
            ax_fowt_disp_pitch = subplot(2,3,5);    
            ax_fowt_disp_yaw = subplot(2,3,6);    
        end
    
        hold(ax_fowt_disp_surge, 'on')
        plot(ax_fowt_disp_surge, data.time, data.surge, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_surge,'Time (s)')
        ylabel(ax_fowt_disp_surge,'Surge (m)')
        set(ax_fowt_disp_surge, 'fontsize', sizeOfFont)
        
        hold(ax_fowt_disp_sway, 'on')
        plot(ax_fowt_disp_sway, data.time, data.sway, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_sway,'Time (s)')
        ylabel(ax_fowt_disp_sway,'Sway (m)')
        set(ax_fowt_disp_sway, 'fontsize', sizeOfFont)

        hold(ax_fowt_disp_heave, 'on')
        plot(ax_fowt_disp_heave, data.time, data.heave, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_heave,'Time (s)')
        ylabel(ax_fowt_disp_heave,'Heave (m)')
        set(ax_fowt_disp_heave, 'fontsize', sizeOfFont)
        
        hold(ax_fowt_disp_roll, 'on')
        plot(ax_fowt_disp_roll, data.roll*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_roll,'Time (s)')
        ylabel(ax_fowt_disp_roll,'Roll (deg)')
        set(ax_fowt_disp_roll, 'fontsize', sizeOfFont)
        
        hold(ax_fowt_disp_pitch, 'on')
        plot(ax_fowt_disp_pitch, data.time, data.pitch*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_pitch,'Time (s)')
        ylabel(ax_fowt_disp_pitch,'Pitch (deg)')
        set(ax_fowt_disp_pitch, 'fontsize', sizeOfFont)
        
        hold(ax_fowt_disp_yaw, 'on')
        plot(ax_fowt_disp_yaw, data.time, data.yaw*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_yaw,'Time (s)')
        ylabel(ax_fowt_disp_yaw,'Yaw (deg)')            
        set(ax_fowt_disp_yaw, 'fontsize', sizeOfFont)
        
        if ii == numel(flNm)
            hl = legend(legSpec, 'position', [0.05, 0.5, 0, 0]);
            set(hl, 'fontsize', 10)                
        end
    end        
end

