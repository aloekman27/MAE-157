clc
clear
close all

load("Pipe_All.txt")

Pipe1 = Pipe_All(1:10,:);
Pipe2 = Pipe_All(11:20,:);
Pipe3 = Pipe_All(21:30,:);
Pipe4 = Pipe_All(31:40,:);

% Time = [Pipe1(:,1) Pipe2(:,1) Pipe3(:,1) Pipe4(:,1)];
% Temp = [Pipe1(:,3) Pipe2(:,3) Pipe3(:,3) Pipe4(:,3)];
% Temp_celcius = Temp - 273.15;
flowRate = [Pipe1(:,4) Pipe2(:,4) Pipe3(:,4) Pipe4(:,4)];
pressureDrop = [Pipe1(:,5) Pipe2(:,5) Pipe3(:,5) Pipe4(:,5)];

D = [8.15 11.8 17.5 9.93] .* 1e-3;
L = [0.781 1.143 1.98 0.337];
K = 5244;



rho = 996; 
viscocity = 0.000798;

vol_fr = ((flowRate./K)./264.172);

area = pi/4.*(D.^2);

velocity = vol_fr ./ area;

% Experimental Friction Factor
fric_factor = (pressureDrop .* (D./L))./(0.5.*rho.*velocity.^2); % y_i

% Experimental Reynolds Number
Re_d = (rho .* velocity .* D)./(viscocity);

%Pipe 1 2 3 Petukhov
f_P = ((0.79.*log(Re_d(:,1:3)))-1.64).^-2; % theoretical y_hat

%Pipe 4 
h = 0.305e-3; %rib height
P = 3.08e-3; %pitch

epsilon = h*exp(3.4 - 0.42*(P/h)^0.46);

% Pipe 4 Nikurdase Equation
f_N = (1.74 + 2.*log10(D(4)./(2.*epsilon))).^-2;
f_4 = f_N*ones(10,1); 

e_over_D = epsilon/D(4);

% Pipe 4 Colebrook - White Equation (Theoretical)
f_CW = (2.*(log10((e_over_D./3.7) - 5.02./Re_d(:,4).*log10(e_over_D./3.7 + 13./Re_d(:,4))))).^-2;


f_theo = [f_P f_CW]; % theoretical friction factor for all pipes

%% Plot of Experimental vs Theoretical
figure
subplot(2,2,1)
scatter(Re_d(:,1),fric_factor(:,1),'*')
hold on
plot(Re_d(:,1),f_theo(:,1))
hold off
title("Friction Factor of Pipe 1")
xlabel("Reynolds Number")
ylabel("friction factor, f")
legend("Experimental","Theoretical")


subplot(2,2,2)
scatter(Re_d(:,2),fric_factor(:,2),'*')
hold on
plot(Re_d(:,2),f_theo(:,2))
hold off
title("Friction Factor of Pipe 2")
xlabel("Reynolds Number")
ylabel("friction factor, f")
legend("Experimental","Theoretical")

subplot(2,2,3)
scatter(Re_d(:,3),fric_factor(:,3),'*')
hold on
plot(Re_d(:,3),f_theo(:,3))
hold off
title("Friction Factor of Pipe 3")
xlabel("Reynolds Number")
ylabel("friction factor, f")
legend("Experimental","Theoretical")

subplot(2,2,4)
scatter(Re_d(:,4),fric_factor(:,4),'*')
hold on
plot(Re_d(:,4),f_theo(:,4))
hold off
title("Friction Factor of Pipe 4")
xlabel("Reynolds Number")
ylabel("friction factor, f")
legend("Experimental","Theoretical")



%% Error Analysis
f_P1 = log(f_P);
f_CW1 = log(f_CW);
Re_1 = log(Re_d);

figure(2)
plot(Re_1(:,1),f_P1(:,1))
hold on
plot(Re_1(:,2),f_P1(:,2))
plot(Re_1(:,3),f_P1(:,3))
plot(Re_1(:,4),f_CW1)


N = 10;

