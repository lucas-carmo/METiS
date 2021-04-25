% Post processing routine for METiS output 

clear all
close all
clc


flNm = { 
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Jappaku\metis\BIC\180DEG-teste\MOD0_BIC20_180DEG_V00_out.txt'
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Jappaku\metis\BIC\180DEG-teste\MOD0_BIC20_180DEG_V00_out_1.txt'

%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Jappaku\metis\WHI\MOD0_WHI01_180DEG_V00_out.txt'
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Jappaku\metis\WHI\MOD0_WHI01_180DEG_V00_out_9.txt'
            
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Misc\Cylinder-2ndOrder_out.txt';
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Misc\Cylinder-2ndOrder_out_4.txt';

%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Misc\Cylinder-JONSWAP_out.txt';
%             'C:\Users\lucas.henrique\Google Drive\Doutorado\1Misc\Cylinder-JONSWAP_out_8.txt';
            
            'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Cilindro\teste_analitico\bic_out.txt';
            'C:\Users\lucas.henrique\Google Drive\Doutorado\1Testes_Cilindro\teste_analitico\bic_out_5.txt';
        };
    
%===== Choose the output
activeDoFs = [1 1 1 1 1 1];

% Use a containers map for easier iteration and identification.
analysisList = containers.Map;

% FOWT position, velocity and acceleration
analysisList('fowt_disp') = 0;
analysisList('fowt_vel') = 0;
analysisList('fowt_acc') = 0;

% Forces
analysisList('hd_inertia_force') = 0;
analysisList('hd_drag_force') = 1;
analysisList('hd_fk_force') = 0;
analysisList('hd_2nd_force_part1') = 0;
analysisList('hd_force') = 1;
analysisList('hs_force') = 0;
analysisList('ad_hub_force') = 0;
analysisList('total_force') = 0;

% Wave outputs
analysisList('wave_elev') = 1;
analysisList('wave_vel') = 0;
analysisList('wave_acc') = 0;
analysisList('wave_acc_2nd') = 0;
analysisList('wave_press') = 0;
analysisList('wave_pres_2nd') = 1;


%===== Plot style
width4Line = 2;
colors4Plot = num2cell(get(groot,'defaultAxesColorOrder'), 2);
sizeOfFont = 12;
style4Plot = {'-', '-', '-'};


%=========================================================================%

% Check if only 0's or 1's were specified in activeDoFs
if numel(activeDoFs(activeDoFs == 1 | activeDoFs == 0)) ~= 6
    error('Value different from 0 or 1 specified in activeDoF.');
end

% Check the number of subplots based on the number of active DoFs          
numSubPlots = sum(activeDoFs);            

k = keys(analysisList);
val = values(analysisList);

% Save the outputs to two three arrays of cells: plotValue, with the time 
% series that will be plotted, plotTitle, with the title of its figure, 
% and plotLabel, with the ylabels of each subplot.
for ii = 1:numel(flNm)    
    data = readOutFl(flNm{ii});
    fieldsOfData = fields(data);
    
    caa = 1; % Count active analysis - counter for the figures where the analysis is active    
    for jj = 1:1:length(analysisList)
        if val{jj} == 0
            continue;
        end                       
        
        title4plot = k(jj);
        label4plot{1,1} = {'surge'; 'sway'; 'heave'; 'roll'; 'pitch'; 'yaw'};
        
        clear y
        % Displacement, forces, etc
        if strcmp(k{jj}, 'fowt_disp')
            y = [data.surge, data.sway, data.heave, data.roll, data.pitch, data.yaw];            
            y(:, activeDoFs == 0) = [];                     
            
        elseif strcmp(k{jj}, 'fowt_vel')
            y = [data.surge_vel, data.sway_vel, data.heave_vel, data.roll_vel, data.pitch_vel, data.yaw_vel];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'fowt_acc')
            y = [data.surge_acc, data.sway_acc, data.heave_acc, data.roll_acc, data.pitch_acc, data.yaw_acc];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hd_inertia_force')
            y = [data.hd_inert_force_1, data.hd_inert_force_2, data.hd_inert_force_3, data.hd_inert_force_4, data.hd_inert_force_5, data.hd_inert_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hd_drag_force')
            y = [data.hd_drag_force_1, data.hd_drag_force_2, data.hd_drag_force_3, data.hd_drag_force_4, data.hd_drag_force_5, data.hd_drag_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hd_fk_force')
            y = [data.hd_fk_force_1, data.hd_fk_force_2, data.hd_fk_force_3, data.hd_fk_force_4, data.hd_fk_force_5, data.hd_fk_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hd_2nd_force_part1')
            y = [data.hd_2nd_force_p1_1, data.hd_2nd_force_p1_2, data.hd_2nd_force_p1_3, data.hd_2nd_force_p1_4, data.hd_2nd_force_p1_5, data.hd_2nd_force_p1_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hd_force')
            y = [data.hd_force_1, data.hd_force_2, data.hd_force_3, data.hd_force_4, data.hd_force_5, data.hd_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'hs_force')
            y = [data.hs_force_1, data.hs_force_2, data.hs_force_3, data.hs_force_4, data.hs_force_5, data.hs_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'ad_hub_force')
            y = [data.ad_hub_force_1, data.ad_hub_force_2, data.ad_hub_force_3, data.ad_hub_force_4, data.ad_hub_force_5, data.ad_hub_force_6];
            y(:, activeDoFs == 0) = [];
            
        elseif strcmp(k{jj}, 'total_force')
            y = [data.total_force_1, data.total_force_2, data.total_force_3, data.total_force_4, data.total_force_5, data.total_force_6];
            y(:, activeDoFs == 0) = [];  
        end 
        
        % The wave outputs are treated a little differently from the ones
        % below, since there can be several wave locations. Besides, some 
        % of them are scalars (elevation and pressure) while others are 
        % vectors (velocity and acceleration).
        if strcmp(k{jj}, 'wave_elev') || strcmp(k{jj}, 'wave_press') || strcmp(k{jj}, 'wave_pres_2nd')
            waveLocation = find(contains(fieldsOfData, k{jj})==1);
            for ww = numel(waveLocation):-1:1
                y(:,1,ww) = data.(fieldsOfData{waveLocation(ww)});
                title4plot{ww} = fieldsOfData(waveLocation(ww));
                label4plot{ww} = {' '; ''; ''};
            end
            
        elseif strcmp(k{jj}, 'wave_vel') || strcmp(k{jj}, 'wave_acc') || strcmp(k{jj}, 'wave_acc_2nd')
            waveLocation_x = find(contains(fieldsOfData, k{jj})==1 & contains(fieldsOfData, '_x')==1); 
            waveLocation_y = find(contains(fieldsOfData, k{jj})==1 & contains(fieldsOfData, '_y')==1); 
            waveLocation_z = find(contains(fieldsOfData, k{jj})==1 & contains(fieldsOfData, '_z')==1); 
            for ww = numel(waveLocation_x):-1:1
                y(:,:,ww) = [data.(fieldsOfData{waveLocation_x(ww)}), data.(fieldsOfData{waveLocation_y(ww)}), data.(fieldsOfData{waveLocation_z(ww)})];
                title4plot{ww} = strrep(fieldsOfData{waveLocation_x(ww)}, '_x', '');
                label4plot{ww} = {'x'; 'y'; 'z'};
            end
        end                
                
        % Only wave quantities may have the third dimension of 'y' larger than 1
        % The others will iterate only once
        for ww = 1:size(y,3)             
            % Create figures and axes
            if ii == 1                           
                figVec(caa) = figure('name', char(title4plot{ww}), 'color', 'w', 'units','normalized','outerposition',[0 0 1 1]);
                title(strrep(title4plot{ww}, '_', '\_'))

                % Check the number of subplots            
                numSubPlots = size(y,2);
                if numSubPlots == 0 && ~contains(k(jj), 'wave')
                    error('At least one DoF should be active when motions or forces are analyzed.');
                end

                % Create distribution of subplots based on the number required
                if numSubPlots >= 4
                    iSubPlot = 2;
                    jSubPlot = 3;

                else
                    iSubPlot = 1;
                    jSubPlot = numSubPlots;
                end

                % Create the axis of eache subplots
                for indSubPlot = 1:numSubPlots
                    axArray{caa, indSubPlot} = subplot(iSubPlot,jSubPlot,indSubPlot);
                end                  
            end                    
                
            % Plot to the corresponding subplot        
            for indSubPlot = 1:size(y,2)
                hold(axArray{caa, indSubPlot}, 'on')
                plot(axArray{caa, indSubPlot}, data.time, y(:,indSubPlot, ww), style4Plot{ii}, 'LineWidth', width4Line, 'color', colors4Plot{ii})
                xlabel(axArray{caa, indSubPlot}, 'Time (s)')
                ylabel(axArray{caa, indSubPlot}, label4plot{ww}{indSubPlot})
                set(axArray{caa, indSubPlot}, 'fontsize', sizeOfFont)                
            end
            caa = caa + 1;
        end        
    end
end


function [f, Yo, Y] = fftamp(x, Fs, graf)
    L = length(x);
    NFFT = 2^nextpow2(L); 
    y = fft(x,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Yo = 2*(y(1:NFFT/2+1));
    Y = 2*abs(y(1:NFFT/2+1));

    % Plot single-sided amplitude spectrum.
    if graf==1
        figure;
        plot(f,2*abs(y(1:NFFT/2+1))); 
        title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
        grid on;
        box on;
    end
end