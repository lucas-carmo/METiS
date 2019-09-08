% Teste pra implementar JONSWAP no METiS

clear all
close all
clc

% Simulation parameters
Tmax = 10800;
dt = 0.25;
t = 0:dt:Tmax;



% JONSWAP parameters
Hs = 2;
Tp = 20;
gamma = 3.3;
dir = 45;
wlow = 0.18;
whigh = 0.92;

% Plot characteristics
width4Line = 2;
size4Font = 25;
weight4Font = 'bold';
aspRatio = [3 1 1];
color4Line = [0, 0.4470, 0.7410];
size4Marker = 10;

addpath('C:\Users\lucas.henrique\Google Drive\Doutorado\1Casos_teste\export-fig')

%========================================================================
dw = 2*pi/Tmax;
wmax = pi/dt;


% Calcula o Jonswap
w = 0.01:dw:wmax;
wp = 2*pi/Tp;

sigma = 0.09*ones(size(w));
sigma(w <= wp) = 0.07;
A = exp( -0.5 * ((w/wp-1)./sigma).^2 );

Sw = 0.3125 * Hs^2 * Tp * (w/wp).^(-5) .* exp(-1.25 * (wp./w).^4 ) .* (1-0.287*log(gamma)) .* gamma.^A;

% figure('color', 'w')
% plot(w,Sw)


% Calcula as amplitudes
amp = 0 * (wlow:dw:whigh);
omega = zeros(size(amp));
Nwaves = length(amp);

ii = 1;
for wloop = wlow : dw : whigh    
    if wloop <= wp
        sigma = 0.07;
    else
        sigma = 0.09;
    end
    A = exp( -0.5 * ((wloop/wp-1)/sigma)^2 );
    S = 0.3125 * Hs^2 * Tp * (wloop/wp)^(-5) * exp(-1.25 * (wp/wloop).^4 ) * (1-0.287*log(gamma)) * gamma^A;
    omega(ii) = wloop;
    amp(ii) = sqrt(2 * S * dw);    
    ii = ii + 1;
end

% Phase is a randon number between -pi and pi
phase = -pi + (pi + pi) * rand(size(amp));


y = zeros(size(t));
for ii = 1:length(amp)
    y = y + amp(ii)*cos(omega(ii)*t + phase(ii));
end


%===
%=== PLOTA FIGURA NO DOMINIO DO TEMPO
%===
outFig = figure('color', 'w', 'units','normalized','outerposition',[0 0 1 1]);
plot(t, y, 'lineWidth', width4Line)
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('First-order wave elevation (m)')
set(gca, 'fontSize', size4Font);
% export_fig(outFig, 'ilust_JONSWAP_timeSeries.pdf');


%===
%=== PLOTA ESPECTRO
%===

% Testar se bate com o espectro de entrada
[Sres, w2] = pwelch(y,[],[],[],1/dt);  
Sres = Sres/(2*pi); % passar o espectro para rad/s 
w2 = 2*pi*w2;


outFig = figure('color', 'w', 'units','normalized','outerposition',[0 0 1 1]);


plot(w2,Sres, 'b', 'lineWidth', width4Line)
hold on
plot(w,Sw, 'r', 'lineWidth', width4Line)
hold on
plot([wlow whigh], [0 0], 'o', 'markerSize', size4Marker, 'Color', 'k', 'MarkerFace', 'k')  
xlim([0 2*whigh])
xlabel('Frequency (rad/s)')
ylabel('Spectral density (m^2.s)')

legend('From time series', 'JONSWAP')
set(gca, 'fontSize', size4Font);

% export_fig(outFig, 'ilust_JONSWAP_espectro.pdf');


