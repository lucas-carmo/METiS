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

flNm = ['..' flSep 'test' flSep 'output_9' flSep 'output.txt'];

data = importdata(flNm);




%%%==
%%%== Teste forca cilindro
%%%==
time = data.data(:,1)';

% Cylinder
n1 = [-5,0,-8.66025 ];
n2 = [5,0,8.66025];
Cm = 1.6;
Cd = 0.5;
D  = 3;
L  = norm(n2-n1);


% Environment
g = 9.81;
h = 2000;
rho = 1025;

%  Wave
A = 1;
omega = 2*pi/10;
direc = 30;

fk = @(x) omega^2/g - x*tanh(x*h);
k = fsolve(fk, 0.1);        

   

% Analitico
H = min(n2(3),n1(3));
alpha = acos( dot( (n2-n1)/L, [0,0,1]) );
a1    = cosd(direc) * cos(alpha);
a2    = sind(direc);
b1    = -sin(alpha);
b2    = 0;

Q = sqrt( ... 
        cos(omega*time).^2 * (a1^2 + a2^2 ) +...
        sin(omega*time).^2 * (b1^2 + b2^2 ) +...
        -2 * cos(omega*time) .* sin(omega*time) * ( a1*b1 + a2*b2 )...
        )';

Fi_analitica = - rho * Cm * (pi/4) * D^2 * ( A/k ) * omega^2 * ( ( 1-exp(k*H) ) / cos(alpha)) * ...
                    [ ( cosd(direc) * cos(alpha) * sin(omega*time)' - sin(alpha) * cos(omega*time)' ) * cos(alpha) ...
                        sind(direc) * sin(omega*time)' ...
                     -( cosd(direc) * cos(alpha) * sin(omega*time)' - sin(alpha) * cos(omega*time)' ) * sin(alpha)
                    ];

Fd_analitica = 0.5 * rho * Cd * D * (A^2/(2*k)) * omega^2 * ( ( 1-exp(2*k*H) )/ cos(alpha) ) * ...
                    repmat( Q , 1 , 3 ) .* ...    
                    [ ( cosd(direc) * cos(alpha) * cos(omega*time)' + sin(alpha) * sin(omega*time)' ) * cos(alpha) ...
                        sind(direc) * cos(omega*time)' ...
                     -( cosd(direc) * cos(alpha) * cos(omega*time)' + sin(alpha) * sin(omega*time)' ) * sin(alpha)
                    ];
                
F_analitica = Fi_analitica + Fd_analitica;


figure
subplot(3,1,1)
plot(time, F_analitica(:,1), 'Color', 'r', 'linewidth', 2)
title(['\beta = ' num2str(direc) 'deg - \alpha = ' num2str(rad2deg(alpha)) 'deg'])
hold on
plot(data.data(:,1),  data.data(:,2), '--', 'Color', 'b', 'linewidth', 2)
hl = legend('Analitico', 'C++');
set(hl, 'FontSize', 20)
set(gca, 'FontSize', 20)
ylabel('Fx (N)');

subplot(3,1,2)
plot(time, F_analitica(:,2), 'Color', 'r', 'linewidth', 2)
hold on
plot(data.data(:,1),  data.data(:,3), '--', 'Color', 'b', 'linewidth', 2)
set(gca, 'FontSize', 20)
ylabel('Fy (N)');

subplot(3,1,3)
plot(time, F_analitica(:,3), 'Color', 'r', 'linewidth', 2)
hold on
plot(data.data(:,1),  data.data(:,4), '--', 'Color', 'b', 'linewidth', 2)
set(gca, 'FontSize', 20)
ylabel('Fz (N)');
xlabel('t (s)')


set(gcf, 'color', 'w');



% % % % %%%==
% % % % %%%== Teste velocidade onda
% % % % %%%==
% % % % t = data.data(:,1);
% % % % 
% % % % g = 9.81;
% % % % h = 200;
% % % % x = 10;
% % % % y = 20;
% % % % z = -0.5;
% % % % 
% % % % H = [2 1 1]';
% % % % omega = [2*pi/10 1 2*pi*10]';
% % % % beta = [45 60 60]';
% % % % 
% % % % k = zeros(length(omega),1);
% % % % khz_xy = zeros(length(omega),1);
% % % % khz_z = zeros(length(omega),1);
% % % % 
% % % % for ii = 1:length(omega)
% % % %     fk = @(x) omega(ii)^2/g - x*tanh(x*h);
% % % %     k(ii) = fsolve(fk, 0.1);        
% % % %     
% % % %     if (k(ii)*h >= 10)
% % % %         khz_xy(ii) = exp(k(ii)*z);
% % % %         khz_z(ii) = khz_xy(ii);
% % % %     else
% % % %         khz_xy(ii) = cosh(k(ii) * (z + h)) / sinh(k(ii)*h);
% % % %         khz_z(ii) = sinh(k(ii) * (z + h)) / sinh(k(ii)*h);
% % % %     end    
% % % %     
% % % % end
% % % % 
% % % % % Essa porra de baixo nao ta funcionando na versao R2014a por causa do cos( omega.*t' )
% % % % % u = sum(omega .* H/2 .* khz_xy .* cos( omega.*t' ) .* cosd(beta), 1);
% % % % % v = sum(omega .* H/2 .* khz_xy .* cos( omega*t' ) .* sind(beta));
% % % % % w = sum(-omega .* H/2 .* khz_z .* sin( omega*t' ));
% % % % u = 0;
% % % % v = 0;
% % % % w = 0;
% % % % a1 = 0;
% % % % a2 = 0;
% % % % a3 = 0;
% % % % for ii = 1:length(omega)
% % % %     u = u + omega(ii) * H(ii)/2 .* khz_xy(ii) * sin( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y ) .* cosd(beta(ii));
% % % %     v = v + omega(ii) * H(ii)/2 .* khz_xy(ii) * sin( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y ) .* sind(beta(ii));
% % % %     w = w + omega(ii) * H(ii)/2 .* khz_z(ii) * cos( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y );
% % % %      
% % % %     a1 = a1 + omega(ii)^2 * H(ii)/2 * khz_xy(ii) * cos( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y )*cosd(beta(ii));
% % % %     a2 = a2 + omega(ii)^2 * H(ii)/2 * khz_xy(ii) * cos( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y )*sind(beta(ii));
% % % %     a3 = a3 - omega(ii)^2 * H(ii)/2 * khz_z(ii) * sin( omega(ii)*t - k(ii)*cosd(beta(ii))*x - k(ii)*sind(beta(ii))*y );
% % % % end
% % % % 
% % % % figure
% % % % subplot(3,1,1)
% % % % plot(t, u, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,2), 'o', 'Color', 'b')
% % % % 
% % % % subplot(3,1,2)
% % % % plot(t, v, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,3), 'o', 'Color', 'b')
% % % % 
% % % % subplot(3,1,3)
% % % % plot(t, w, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,4), 'o', 'Color', 'b')
% % % % 
% % % % 
% % % % 
% % % % figure
% % % % subplot(3,1,1)
% % % % plot(t, a1, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,5), 'o', 'Color', 'b')
% % % % 
% % % % subplot(3,1,2)
% % % % plot(t, a2, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,6), 'o', 'Color', 'b')
% % % % 
% % % % subplot(3,1,3)
% % % % plot(t, a3, 'Color', 'r', 'linewidth', 2)
% % % % hold on
% % % % plot(data.data(1:5:end,1),  data.data(1:5:end,7), 'o', 'Color', 'b')


