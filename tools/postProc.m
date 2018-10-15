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

flNm = ['..' flSep 'test' flSep 'output_13' flSep 'output.txt'];

data = importdata(flNm);


% Teste velocidade onda
t = data.data(:,1);

g = 9.81;
h = 200;
z = 0;

H = [2 4 4]';
omega = [2*pi/10 1 2*pi*10]';
beta = [45 60 60]';

k = zeros(length(omega),1);
khz_xy = zeros(length(omega),1);
khz_z = zeros(length(omega),1);

for ii = 1:length(omega)
    fk = @(x) omega(ii)^2/g - x*tanh(x*h);
    k(ii) = fsolve(fk, 0.1);        
    
    if (k(ii)*h >= 10)
        khz_xy(ii) = exp(k(ii)*z);
        khz_z(ii) = khz_xy(ii);
    else
        khz_xy(ii) = cosh(k(ii) * (z + h)) / sinh(k(ii)*h);
        khz_z(ii) = sinh(k(ii) * (z + h)) / sinh(k(ii)*h);
    end    
    
end

% Essa porra de baixo nao ta funcionando na versao R2014a por causa do cos( omega.*t' )
% u = sum(omega .* H/2 .* khz_xy .* cos( omega.*t' ) .* cosd(beta), 1);
% v = sum(omega .* H/2 .* khz_xy .* cos( omega*t' ) .* sind(beta));
% w = sum(-omega .* H/2 .* khz_z .* sin( omega*t' ));
u = 0;
v = 0;
w = 0;
for ii = 1:length(omega)
    u = u + omega(ii) * H(ii)/2 .* khz_xy(ii) * cos( omega(ii)*t ) .* cosd(beta(ii));
    v = v + omega(ii) * H(ii)/2 .* khz_xy(ii) * cos( omega(ii)*t ) .* sind(beta(ii));
    w = w - omega(ii) * H(ii)/2 .* khz_z(ii) * sin( omega(ii)*t );
end

figure

subplot(3,1,1)
plot(t, u, 'Color', 'r', 'linewidth', 2)
hold on
plot(data.data(1:5:end,1),  data.data(1:5:end,2), 'o', 'Color', 'b')

subplot(3,1,2)
plot(t, v, 'Color', 'r', 'linewidth', 2)
hold on
plot(data.data(1:5:end,1),  data.data(1:5:end,3), 'o', 'Color', 'b')

subplot(3,1,3)
plot(t, w, 'Color', 'r', 'linewidth', 2)
hold on
plot(data.data(1:5:end,1),  data.data(1:5:end,4), 'o', 'Color', 'b')





