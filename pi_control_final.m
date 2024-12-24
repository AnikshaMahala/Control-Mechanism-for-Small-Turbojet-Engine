clc; clear; 
format short e 

% Parameters and initialization
i=1; m=1; n=1; c=1; q=1; kl=0;
Pt_2=93380;  % Pascal 
Tt_2(1)=303;    % Kelvin 
Pt_ref=100000; % Pascal (reference) 
Tt_ref=298;  % Kelvin (reference) 
sigma_comb=0.95;  % Combustion chamber pressure ratio
nu_comb=0.6;  % Combustion chamber efficiency
Pi_n=0.95;  % Nozzle pressure ratio
M9=1;  % Choke condition
M0=0;  % Assumption
H=0;  % Height

% Additional initialization
teta=(Tt_2(1)/Tt_ref);  % Dimensionless temperature parameter
delta=(Pt_2/Pt_ref); % Dimensionless pressure parameter
N(1)=65000;  % RPM   
I=1.787e-4;  % Moment of Inertia _from CATIA [kg*m^2] 
eff_t=0.85;  % Turbine efficiency
gama_c=1.4;  % Compressor specific heat ratio
gama_t=1.35; % Turbine specific heat ratio
Qf=43250e3;  % [J/kg] Lower calorific value of fuel
R=286.9;  % [J/(kg*K)] Individual gas constant
V_comb=6.687233e-4;  % Volume of the combustion chamber

% Compressor-related parameters
PI_C=9.06e-21*N(1)^4 - 1.365e-15*N(1)^3 + 3.3806e-10*N(1)^2 - 9.7556e-6*N(1) + 1.1929;
Cp_a(1)=-1.9747e-028*N(1)^6 + 8.3860e-023*N(1)^5 - 1.4447e-017*N(1)^4 + 1.2943e-012*N(1)^3 - 6.2674e-008*N(1)^2 + 1.5743e-003*N(1) + 9.8846e+002;    
Tt_4(1)=5.5099e-017*N(1)^4 - 1.1969e-011*N(1)^3 + 9.7144e-007*N(1)^2 - 3.7423e-002*N(1) + 1.4273e+003;

% Initialize mass flow rates
mf_dot(1)=0.00456;
Kp=0.00000456/10000; 
Ki=0.00000000035;

MODD=2; % MODD value based on your code
A_0=0.001378; 
A_9=0.00302; 

% Delta "t"
h=0.001;    
E(1)=0; 

%-----------BETA LINE METHOD---------------------------------------- 
for j=1:(1.5/h) 
    N_Ref(j)=85000; 
    N_R(j)=N(j)/(teta)^0.5;  % Corrected RPM
    
    % Ensure correct initialization of Tt_2 and Pt_2
    Tt_2(j+1)=Tt_2(j); 
    Pt_2(j+1)=Pt_2(j); 
    
    % Map data initialization (simplified for readability)
    PI_C_m=[0 1.560 1.648 1.678 0 0 0 0 0 0 0 0  
            0 2.019 2.172 2.282 2.337 0 0 0 0 0 0 0   
            0 0 2.545 2.750  2.938 3.067 3.121 0 0 0 0 0   
            0  0 0 0 3.219 3.438 3.676 3.843 3.955 0 0 0  
            0 0 0 0 0 0 3.861 4.090 4.317 4.542 4.695 0]; 
    N_cor=[61000 81000 97000 111000 123000 ]; 
    
    % Compressor interpolation using Beta Line Method
    while m==1 
        if N_cor(i) <= N_R(j) && N_R(j) <= N_cor(i+1) 
            k=(N_cor(i+1)-N_R(j))/(N_cor(i+1)-N_cor(i)); 
            m=0; 
            t=i; 
        end
        i=i+1; 
    end
    i=1; 
    
    % Interpolation of compressor map values
    for q = 1 : 12 
        N_R_matrix(q) = PI_C_m(t+1,q) - (PI_C_m(t+1,q) - PI_C_m(t,q)) * k; 
    end
    
    % Calculation of Beta for mass flow rate calculations
    while n==1 
        if N_R_matrix(c) <= PI_C(j) && PI_C(j) <= N_R_matrix(c+1) 
            Beta = c - (-1 * (N_R_matrix(c) - PI_C(j))) / (N_R_matrix(c) - N_R_matrix(c+1)); 
            n=0; 
        end
        c=c+1; 
    end
    c=1; 
    
    % Mass Flow Rate Calculations
    x = [1 2 3 4 5 6 7 8 9 10 11 12]; 
    y = [61000 81000 97000 111000 123000 ]; 
    z = [0 0.2777 0.2185 0.1387 0 0 0 0 0 0 0 0  
         0 0.4155 0.3758 0.3226 0.2531 0 0 0 0 0 0 0  
         0 0 0.4904 0.4638 0.4358 0.3894 0.318 0 0 0 0 0 
         0 0 0 0 0.5206 0.5045 0.4881 0.4514 0.3988 0 0 0  
         0 0 0 0 0 0 0.5435 0.5271 0.5094 0.4918 0.4517 0];
     
    % Efficiency Calculations
    x_1 = [1 2 3 4 5 6 7 8 9 10 11 12];
    y_1 = [61000 81000 97000 111000 123000 ];
    z_1 = [0 0.711 0.759 0.69 0 0 0 0 0 0 0 0  
           0 0.669 0.733 0.76 0.742 0 0 0 0 0 0 0  
           0 0 0.646 0.698 0.735 0.752 0.726 0 0 0 0 0  
           0 0 0 0 0.65 0.68 0.704 0.715 0.706 0 0 0  
           0 0 0 0 0 0 0.62 0.64 0.655 0.664 0.66 0]; 
    
    % Calculate mass flow rate and efficiency
    m_dot_MAP(j)=interp2(x,y,z,Beta,N_R(j)); 
    eff_c(j)=interp2(x_1,y_1,z_1,Beta,N_R(j)); 
    
    %------------------------Compressor--------------------------------- 
    mc_dot(j)=m_dot_MAP(j)*(delta)/((teta)^0.5); 
    Tt_3(j)=Tt_2(j)*(1+((1/eff_c(j))*((PI_C(j)^((gama_c-1)/gama_c))-1))); 
    Pt_3(j)=PI_C(j)*Pt_2(j); 
    Ta(j)=(Tt_2(j)+Tt_3(j))/2; 
    Cp_a(j+1)=1.0189e3 - 0.13784*Ta(j) + 1.9843e-4*Ta(j)^2 + 4.2399e-7*Ta(j)^3 - 3.7632e-10*Ta(j)^4;
    
    %------------------------Combustion Chamber------------------------- 
    mf_dot(j+1)=mf_dot(j); 
    kont=rem(j,MODD); 
    if kont==0 
        E(j)=N_Ref(j)-N(j); 
        mf_dot(j)=mf_dot(j-1)+ (Kp+Ki*MODD).*E(j)+Kp*E(j-1);
        mf_dot(j+1)=mf_dot(j); 
    end 
    
    Tg(j)=(Tt_4(j)+Tt_3(j))/2; 
    Tm(j)=Tg(j); 
    f(j)=mf_dot(j)/mc_dot(j); 
    
    % Beta coefficients for gas
    Bt(j)=-3.59494e2 + 4.5164*Tg(j) + 2.8116e-3*Tg(j)^2 - 2.1709e-5*Tg(j)^3 + 2.8689e-8*Tg(j)^4 - 1.2263e-11*Tg(j)^5; 
    Cp_g(j)=Cp_a(j)+(f(j)/(1+f(j)))*Bt(j); 
    Cv_comb(j)=Cp_g(j)-R;  
    
    Pt_4(j)=sigma_comb*Pt_3(j); 
    Pm(j)=(Pt_4(j)+Pt_3(j))/2; 
    Rho(j)=Pm(j)/(R*Tm(j)); 
    M_comb(j)=Rho(j)*V_comb; 
    
    %------------------Runge Kutta fourth order for Tt_4 ---------------           
    F(j)=(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*Tt_4(j)+Qf*nu_comb*mf_dot(j))/(Cv_comb(j)*M_comb(j)); 
    K1(j)=F(j); 
    K2(j)=(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*(Tt_4(j)+(h/2)*K1(j))+Qf*nu_comb*mf_dot(j))/(Cv_comb(j)*M_comb(j)); 
    K3(j)=(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*(Tt_4(j)+(h/2)*K2(j))+Qf*nu_comb*mf_dot(j))/(Cv_comb(j)*M_comb(j)); 
    K4(j)=(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*(Tt_4(j)+h*K3(j))+Qf*nu_comb*mf_dot(j))/(Cv_comb(j)*M_comb(j)); 
    
    Tt_4(j+1)=Tt_4(j)+(h/6)*(K1(j)+2*K2(j)+2*K3(j)+K4(j)); 
    
    %------------------Runge Kutta fourth order for Pt_4----------------          
    H(j)=Pt_4(j)*(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*Tt_4(j)+Qf*nu_comb*mf_dot(j))/(Tt_4(j)*Cv_comb(j)*M_comb(j)); 
    M1(j)=H(j); 
    M2(j)=(Pt_4(j)+h/2*M1(j))*(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*Tt_4(j)+Qf*nu_comb*mf_dot(j))/(Tt_4(j)*Cv_comb(j)*M_comb(j)); 
    M3(j)=(Pt_4(j)+h/2*M2(j))*(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*Tt_4(j)+Qf*nu_comb*mf_dot(j))/(Tt_4(j)*Cv_comb(j)*M_comb(j)); 
    M4(j)=(Pt_4(j)+h*M3(j))*(mc_dot(j)*Cp_a(j)*Tt_3(j)-(mc_dot(j)+mf_dot(j))*Cp_g(j)*Tt_4(j)+Qf*nu_comb*mf_dot(j))/(Tt_4(j)*Cv_comb(j)*M_comb(j)); 
    
    Pt_4(j+1)=Pt_4(j)+h/6*(M1(j)+2*M2(j)+2*M3(j)+M4(j));   
    Pt_3(j+1)=Pt_4(j+1)/sigma_comb; 
    PI_C(j+1)=Pt_3(j+1)/Pt_2(j); 
    
    %--------------------TURBINE---------------------------------------- 
    PI_t(j)=-4.5415e-029*N(j)^6+1.7103e-023*N(j)^5-2.6255e-018*N(j)^4+2.1127e-013*N(j)^3-9.4373e-009*N(j)^2+2.1571e-004*N(j)-1.0868; 
    Tao_c(j)=Tt_3(j)/Tt_2(j); 
    Teta_2(j)=1; 
    Teta_3(j)=Tt_4(j)/Tt_2(j); 
    Tao_t(j)=1-((Cp_a(j)*Teta_2(j))/(Cp_g(j)*(1+f(j))*Teta_3(j)))*(Tao_c(j)-1); 
    Tt_5(j)=685; 
    
    %--------------------Work Balance----------------------------------- 
    Pt(j)=((mc_dot(j)+ mf_dot(j))*(Tt_4(j)-Tt_5(j))*Cp_g(j)); 
    Pc(j)=mc_dot(j)*Cp_a(j)*(Tt_3(j)-Tt_2(j)); 
    Pf(j)=500; 
    if Pf(j)<0 
        Pf(j)=0; 
    end 
    W(j)=N(j)*pi/30; 
    W(j+1)=sqrt((2*(Pt(j)-Pc(j)-Pf(j))/I)*h+W(j)^2); 
    N(j+1)=W(j+1)*30/pi; 
    
    %------------------THRUST------------------------------------------
    Rho_0(j)=-1.4267e-029*N(j)^6 + 5.8500e-024*N(j)^5 - 9.7835e-019*N(j)^4 + 8.5570e-014*N(j)^3 - 4.1601e-009*N(j)^2 + 1.0194e-004*N(j) + 1.9642e-001; 
    m9_dot(j)=mc_dot(j)+ mf_dot(j); 
    T_9(j)=8.6554e-026*N(j)^6 - 3.4961e-020*N(j)^5 + 5.8197e-015*N(j)^4 - 5.0939e-010*N(j)^3 + 2.4663e-005*N(j)^2 - 6.2630e-001*N(j) + 7.2095e+003; 
    V9(j)=(gama_t*R*T_9(j))^0.5; 
    V0(j)=mc_dot(j)/(Rho_0(j)*A_0); 
    F(j)=m9_dot(j)*V9(j)-mc_dot(j)*V0(j); 
end 

% Plot results
subplot(1,3,1); plot(N, 'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('N (rpm)');
subplot(1,3,2); plot(Tt_4, 'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('Tt_4 (K)')
subplot(1,3,3); plot(F, 'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('Thrust(N)')








