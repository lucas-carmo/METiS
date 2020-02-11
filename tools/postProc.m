% Post processing routine for METiS output 

clear all
close all
clc


flNm = { '/home/lucas/Nextcloud/Doutorado/Casos teste/preliminares/old/Verify Aero/verifyAero_caso3_out.txt'
         '/home/lucas/Work/METiS/test/verifyAero/verifyAero_case3_out_2.txt'
        };

legSpec = { 'Velho'
            'Novo'
          };
    
%===== Choose the output
activeDoFs = [1 1 1 1 1 1];
fowt_disp = 0;
fowt_vel = 0;
fowt_acc = 0;

hd_inertia_force = 0;
hd_drag_force = 0;
hd_fk_force = 0;
hd_force = 0;
hs_force = 0;
ad_hub_force = 1;
total_force = 0;

% Wave outputs
wave_elev = 0;
wave_vel = 0;
wave_acc = 0;
wave_press = 0;



%===== Plot style
width4Line = 2;
colors4Plot = num2cell(get(groot,'defaultAxesColorOrder'), 2);
sizeOfFont = 12;


%=========================================================================%

% Check if only 0's or 1's were specified in activeDoFs
if numel(activeDoFs(activeDoFs == 1 | activeDoFs == 0)) ~= 6
    error('Value different from 0 or 1 specified in activeDoF.');
end

for ii = 1:numel(flNm)
    data = readOutFl(flNm{ii});
    
    if ad_hub_force == 1
        if ii == 1
           fig_ad_hub_force_forces = figure('name', 'Hub aero forces' ,'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);
           ax_fowt_disp_fx = subplot(3,1,1);
           ax_fowt_disp_fy = subplot(3,1,2);
           ax_fowt_disp_fz = subplot(3,1,3);
           
           fig_ad_hub_force_moments = figure('name', 'Hub aero moments' ,'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);
           ax_fowt_disp_mx = subplot(3,1,1);
           ax_fowt_disp_my = subplot(3,1,2);
           ax_fowt_disp_mz = subplot(3,1,3);                     
        end        
        
        hold(ax_fowt_disp_fx, 'on')
        plot(ax_fowt_disp_fx, data.time, data.ad_hub_force_1, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_fx, 'Time (s)')
        ylabel(ax_fowt_disp_fx, 'Fx (N)')
        title(ax_fowt_disp_fx, 'Hub aero force')  
        fprintf('Fx media (%d): %.4e\n', ii, mean(data.ad_hub_force_1));
        fprintf('Fx max (%d): %.4e\n', ii, max(data.ad_hub_force_1));
        fprintf('Fx min (%d): %.4e\n', ii, min(data.ad_hub_force_1));
        fprintf('Fx amp (%d): %.4e\n', ii, max(data.ad_hub_force_1)-min(data.ad_hub_force_1));
        
        hold(ax_fowt_disp_fy, 'on')
        plot(ax_fowt_disp_fy, data.time, data.ad_hub_force_2, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_fy, 'Time (s)')
        ylabel(ax_fowt_disp_fy, 'Fy (N)')

        hold(ax_fowt_disp_fz, 'on')
        plot(ax_fowt_disp_fz, data.time, data.ad_hub_force_3, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_fz, 'Time (s)')
        ylabel(ax_fowt_disp_fz, 'Fz (N)')
                 
        hold(ax_fowt_disp_mx, 'on')
        plot(ax_fowt_disp_mx, data.time, data.ad_hub_force_4, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_mx, 'Time (s)')
        ylabel(ax_fowt_disp_mx, 'Mx (N)')
        title(ax_fowt_disp_mx, 'Hub aero moments')
        
        hold(ax_fowt_disp_my, 'on')
        plot(ax_fowt_disp_my, data.time, data.ad_hub_force_5, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_my, 'Time (s)')
        ylabel(ax_fowt_disp_my, 'My (N)')

        hold(ax_fowt_disp_mz, 'on')
        plot(ax_fowt_disp_mz, data.time, data.ad_hub_force_6, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_mz, 'Time (s)')
        ylabel(ax_fowt_disp_mz, 'Mz (N)')
        
        if ii == numel(flNm)
            hl = legend(legSpec, 'position', [0.05, 0.5, 0, 0]);
            set(hl, 'fontsize', 10)                
        end
    end
    
    
    if fowt_disp == 1
        if ii == 1
           fig_fowt_disp = figure('name', 'FOWT displacement', 'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);
            
            % Check the number of subplots            
            numSubPlots = sum(activeDoFs);            
            if numSubPlots == 0
                error('At least one DoF should be active when fowt_disp==1');
            end
            
            % Create distribution of subplots based on the number required
            if numSubPlots >= 4
                iSubPlot = 2;
                jSubPlot = 3;
                
            else
                iSubPlot = 1;
                jSubPlot = numSubPlots;
            end
            
            % Create the axis of each subplot            
            iSubPlots = 1;
            if activeDoFs(1) == 1
                ax_fowt_disp_surge = subplot(iSubPlot,jSubPlot,iSubPlots);  
                iSubPlots = iSubPlots+1;
            end

            if activeDoFs(2) == 1
                ax_fowt_disp_sway = subplot(iSubPlot,jSubPlot,iSubPlots);    
                iSubPlots = iSubPlots+1;
            end

            if activeDoFs(3) == 1
                ax_fowt_disp_heave = subplot(iSubPlot,jSubPlot,iSubPlots);    
                iSubPlots = iSubPlots+1;
            end

            if activeDoFs(4) == 1
                ax_fowt_disp_roll = subplot(iSubPlot,jSubPlot,iSubPlots);    
                iSubPlots = iSubPlots+1;
            end

            if activeDoFs(5) == 1
                ax_fowt_disp_pitch = subplot(iSubPlot,jSubPlot,iSubPlots);    
                iSubPlots = iSubPlots+1;
            end

            if activeDoFs(6) == 1
                ax_fowt_disp_yaw = subplot(iSubPlot,jSubPlot,iSubPlots);                                            
                iSubPlots = iSubPlots+1;
            end            
        end
            
        if activeDoFs(1) == 1
            hold(ax_fowt_disp_surge, 'on')
            plot(ax_fowt_disp_surge, data.time, data.surge, 'LineWidth', width4Line, 'color', colors4Plot{ii})
            xlabel(ax_fowt_disp_surge,'Time (s)')
            ylabel(ax_fowt_disp_surge,'Surge (m)')
            set(ax_fowt_disp_surge, 'fontsize', sizeOfFont)
        end
        
        if activeDoFs(2) == 1
            hold(ax_fowt_disp_sway, 'on')
            plot(ax_fowt_disp_sway, data.time, data.sway, 'LineWidth', width4Line, 'color', colors4Plot{ii})
            xlabel(ax_fowt_disp_sway,'Time (s)')
            ylabel(ax_fowt_disp_sway,'Sway (m)')
            set(ax_fowt_disp_sway, 'fontsize', sizeOfFont)
        end

        if activeDoFs(3) == 1
            hold(ax_fowt_disp_heave, 'on')
            plot(ax_fowt_disp_heave, data.time, data.heave, 'LineWidth', width4Line, 'color', colors4Plot{ii})
            xlabel(ax_fowt_disp_heave,'Time (s)')
            ylabel(ax_fowt_disp_heave,'Heave (m)')
            set(ax_fowt_disp_heave, 'fontsize', sizeOfFont)
        end
        
        if activeDoFs(4) == 1
        hold(ax_fowt_disp_roll, 'on')
        plot(ax_fowt_disp_roll, data.roll*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
        xlabel(ax_fowt_disp_roll,'Time (s)')
        ylabel(ax_fowt_disp_roll,'Roll (deg)')
        set(ax_fowt_disp_roll, 'fontsize', sizeOfFont)
        end
        
        if activeDoFs(5) == 1
            hold(ax_fowt_disp_pitch, 'on')
            plot(ax_fowt_disp_pitch, data.time, data.pitch*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
            xlabel(ax_fowt_disp_pitch,'Time (s)')
            ylabel(ax_fowt_disp_pitch,'Pitch (deg)')
            set(ax_fowt_disp_pitch, 'fontsize', sizeOfFont)
        end
        
        if activeDoFs(6) == 1
            hold(ax_fowt_disp_yaw, 'on')
            plot(ax_fowt_disp_yaw, data.time, data.yaw*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
            xlabel(ax_fowt_disp_yaw,'Time (s)')
            ylabel(ax_fowt_disp_yaw,'Yaw (deg)')            
            set(ax_fowt_disp_yaw, 'fontsize', sizeOfFont)
        end
        
        if ii == numel(flNm)
            hl = legend(legSpec, 'position', [0.05, 0.5, 0, 0]);
            set(hl, 'fontsize', 10)                
        end
    end        
end

% % % % Input:
% % % % - Matrices xInput e yInput with size lengthSignal x numPlots. Each column
% % % % is plotted to a different subplot. xInput(:,ii) with yInput(:,ii) in one,
% % % % xInput(:,ii+1) with yInput(:,ii+1) in another, etc.
% % % % - Name for the figure
% % % function plotDoFs(xInput, yInput, figName)
% % %            hFig = figure('name', figName ,'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);                       
% % %             
% % %             % Check the number of subplots            
% % %             if (size(xInput,2) ~= size(yInput,2))
% % %             
% % %             numSubPlots = size(xInput,2);
% % %             
% % %             
% % %             
% % %             if numSubPlots == 0
% % %                 return;
% % %             end
% % %             
% % %             % Create distribution of subplots based on the number required
% % %             if numSubPlots >= 4
% % %                 iSubPlot = 2;
% % %                 jSubPlot = 3;
% % %                 
% % %             else
% % %                 iSubPlot = 1;
% % %                 jSubPlot = numSubPlots;
% % %             end
% % %             
% % %             % Create the axis of each subplot            
% % %             iSubPlots = 1;
% % %             if activeDoFs(1) == 1
% % %                 ax_fowt_disp_surge = subplot(iSubPlot,jSubPlot,iSubPlots);  
% % %                 iSubPlots = iSubPlots+1;
% % %             end
% % % 
% % %             if activeDoFs(2) == 1
% % %                 ax_fowt_disp_sway = subplot(iSubPlot,jSubPlot,iSubPlots);    
% % %                 iSubPlots = iSubPlots+1;
% % %             end
% % % 
% % %             if activeDoFs(3) == 1
% % %                 ax_fowt_disp_heave = subplot(iSubPlot,jSubPlot,iSubPlots);    
% % %                 iSubPlots = iSubPlots+1;
% % %             end
% % % 
% % %             if activeDoFs(4) == 1
% % %                 ax_fowt_disp_roll = subplot(iSubPlot,jSubPlot,iSubPlots);    
% % %                 iSubPlots = iSubPlots+1;
% % %             end
% % % 
% % %             if activeDoFs(5) == 1
% % %                 ax_fowt_disp_pitch = subplot(iSubPlot,jSubPlot,iSubPlots);    
% % %                 iSubPlots = iSubPlots+1;
% % %             end
% % % 
% % %             if activeDoFs(6) == 1
% % %                 ax_fowt_disp_yaw = subplot(iSubPlot,jSubPlot,iSubPlots);                                            
% % %                 iSubPlots = iSubPlots+1;
% % %             end            
% % %         end
% % %             
% % %         if activeDoFs(1) == 1
% % %             hold(ax_fowt_disp_surge, 'on')
% % %             plot(ax_fowt_disp_surge, data.time, data.surge, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %             xlabel(ax_fowt_disp_surge,'Time (s)')
% % %             ylabel(ax_fowt_disp_surge,'Surge (m)')
% % %             set(ax_fowt_disp_surge, 'fontsize', sizeOfFont)
% % %         end
% % %         
% % %         if activeDoFs(2) == 1
% % %             hold(ax_fowt_disp_sway, 'on')
% % %             plot(ax_fowt_disp_sway, data.time, data.sway, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %             xlabel(ax_fowt_disp_sway,'Time (s)')
% % %             ylabel(ax_fowt_disp_sway,'Sway (m)')
% % %             set(ax_fowt_disp_sway, 'fontsize', sizeOfFont)
% % %         end
% % % 
% % %         if activeDoFs(3) == 1
% % %             hold(ax_fowt_disp_heave, 'on')
% % %             plot(ax_fowt_disp_heave, data.time, data.heave, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %             xlabel(ax_fowt_disp_heave,'Time (s)')
% % %             ylabel(ax_fowt_disp_heave,'Heave (m)')
% % %             set(ax_fowt_disp_heave, 'fontsize', sizeOfFont)
% % %         end
% % %         
% % %         if activeDoFs(4) == 1
% % %         hold(ax_fowt_disp_roll, 'on')
% % %         plot(ax_fowt_disp_roll, data.roll*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %         xlabel(ax_fowt_disp_roll,'Time (s)')
% % %         ylabel(ax_fowt_disp_roll,'Roll (deg)')
% % %         set(ax_fowt_disp_roll, 'fontsize', sizeOfFont)
% % %         end
% % %         
% % %         if activeDoFs(5) == 1
% % %             hold(ax_fowt_disp_pitch, 'on')
% % %             plot(ax_fowt_disp_pitch, data.time, data.pitch*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %             xlabel(ax_fowt_disp_pitch,'Time (s)')
% % %             ylabel(ax_fowt_disp_pitch,'Pitch (deg)')
% % %             set(ax_fowt_disp_pitch, 'fontsize', sizeOfFont)
% % %         end
% % %         
% % %         if activeDoFs(6) == 1
% % %             hold(ax_fowt_disp_yaw, 'on')
% % %             plot(ax_fowt_disp_yaw, data.time, data.yaw*180/pi, 'LineWidth', width4Line, 'color', colors4Plot{ii})
% % %             xlabel(ax_fowt_disp_yaw,'Time (s)')
% % %             ylabel(ax_fowt_disp_yaw,'Yaw (deg)')            
% % %             set(ax_fowt_disp_yaw, 'fontsize', sizeOfFont)
% % %         end
% % %         
% % %         if ii == numel(flNm)
% % %             hl = legend(legSpec, 'position', [0.05, 0.5, 0, 0]);
% % %             set(hl, 'fontsize', 10)                
% % %         end   
% % % end