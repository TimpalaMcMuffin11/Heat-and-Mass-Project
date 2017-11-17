# Heat-and-Mass-Project
Optimizing Pin Heat Transfer for IC Chip

%% Heat and Mass Project

clear all
close all
clc

% Problem Constants

L_chip = 0.03; % Chip Length (m)
W_chip = 0.03; % Chip Width (m)
A_chip = L_chip*W_chip; % Chip Area (m^2)
L = 0.005; % Board Thickness (m)
k_board = 15; % Board Thermal Conductivity (W/mK)

h = 1000:2000:5000; % Convection Heat Transfer Coefficient (W/m^2K)

% Heat Transfer without Fins

T_s = 80; % Surface Temperature (C)
T_inf = 25:35; % Air Temperature (C)
theta_b = T_s-T_inf; % Differential Temperature (C)
R_tc = 1e-04; % Thermal Contact Resistance (m^2)(K)/W
R1=(1./(h*A_chip));
R2=(R_tc + (L/(k_board*A_chip))+(1./(h*A_chip)));
R_total = 1./((1./R1)+(1./(R2))); % Overall Board Resistance
[X,Y] = meshgrid(theta_b,R_total);

q_max = X./Y; % Heat Transfer (W)

%Table = array2table([T_inf',q_max((h==10),:)',(q_max((h==5000),:)')],'VariableNames',{'T_inf','Convection_10','Convection_5000'});

plot(h,q_max)
title('Heat Transfer vs. Convection Coefficient (without fins)')
xlabel('Convection Coefficient \it h')
ylabel('Heat Transfer \it W')
legend(((num2str(T_inf','T_i_n_f=%d'))),'Location','southeast')

% Design Constants

k_Al = 279; % Thermal Conductivity of Aluminum (W/mK)
k_Cu = 401; % Thermal Conductivity of Copper (W/mK)
L_b = 0.003; % Heat Sink Base Height (m)

% Design Parameters (Variable)

D_p = 0.001:0.001:0.003; % Fin Diameter (m)
L_p = 0.01:0.01:0.03; % Fin Length (m)
N = 5:1:15; % Array

% Geometry Calculations (Dependent on Parameter Variation)

P = pi*D_p; % Fin Perimeter (m)  
A_c = pi/4*(D_p.^2); % Fin Cross Sectional Area (m^2)

L_g = zeros(length(D_p),length(N)); % Length Pin gap (m)
for i = 1:length(N)
  for  j = 1:length(D_p)
    L_g(j,i) = (L_chip-((N(i))*D_p(j)))/(N(i)+1); 
    if(L_g(j,i)<0)
        L_g(j,i)=0;
    end  
  end
end

L_c=zeros(length(L_p),length(D_p)); % Corrected Length (m)
for i = 1:length(D_p)
    for j = 1:length(L_p)
        L_c(j,i) = L_p(j) + (D_p(i)/4); 
    end
end

% Pin Efficiency Calculations

m=zeros(length(h),length(D_p));
for i = 1:length(D_p)
    for j = 1:length(h)
        m(j,i)=sqrt((h(j)*P(i))/(k_Al*A_c(i)));
    end
end

A_f=zeros(numel(L_c),length(D_p)); % Fin Area (m^2)
for i = 1:length(D_p)
    for j = 1:numel(L_c)
        A_f(j,i) = pi*D_p(i)*L_c(j); 
    end
end

A_b = zeros(length(N),length(A_c));
for i = 1:length(A_c)
  for  j = 1:length(N)
    A_b(j,i) = A_chip-(N(j)^2*A_c(i)); % Exposed Base Area (m^2)
  end
end

A_t = zeros(length(N),numel(A_f),numel(A_b));
for i = 1:length(N)
    for j = 1:numel(A_f)
        for k = 1:numel(A_b)
            A_t(j,k,i) = (N(i)*A_f(j))+A_b(k); % Total Surface Area (m^2)
        end
    end
end

n_f = zeros(numel(m),length(L_c));
for i = 1:numel(m)
    for j = 1:length(L_c)
        n_f(j,i) = tanh(m(j)*L_c(i))/(m(j)*L_c(i)); % Fin Efficiency
    end
end

n_o = zeros(length(N),numel(A_f),numel(A_t),numel(n_f));
for i = 1:length(N)
    for j = 1:numel(A_f)
        for k = 1:numel(A_t)
            for u = 1:numel(n_f)
                n_o = 1-((N(i)*A_f(j)/A_t(k))*(1-n_f(u))); % Overall Surface Efficiency
            end
        end
    end
end











% R_to = 1/(n_o*h*A_t); % Fin Thermal Resistance
% q_f = theta_b/R_to; % Fin Heat Transfer (W)
% q_c = q_board+q_f % Total Heat Transfer with Fins (W)


