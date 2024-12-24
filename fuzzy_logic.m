%{
clc;clear; 
format short e 
i=1;m=1;nn=1;c=1;q=1;v=1;kl=0; 
Pt_2=93380;          %pascal 
Tt_2=303;            %Kelvin 
Pt_ref=100000;       %pascal  (from compressor map) 
Tt_ref=298;          %Kelvin  (from compressor map) 
sigma_comb= 0.95;    %Combustion chamber pressor ratio 
nu_comb=0.6;         %Combustion chamber efficiency
Pi_n=0.95;           %Nozzle pressure ratio 
M9=1;                %Choke condition 
M0=0;                %Assumption
H=0;                 %Hight
teta=(Tt_2/Tt_ref);  %Dimensionless temperature parameter 
delta=(Pt_2/Pt_ref); %Dimensionless pressure parameter 
N=65000  ;           %RPM   
I=1.787e-4;          % Moment of Inertia _from CATIA [kg*m^2] 
eff_t=0.85;          %efficiency of the turbine (Can Dedekargınoğlu) 
gama_c=1.4;         %specific heat ratio of compressor 
gama_t=1.35;       %specific heat ratio of turbine   
Qf=43250e3;        % [J/kg]  Lower calorific value of fuel  
R=286.9;            % [J/(kg*K)] Individual gas constant 
V_comb=6.687233e-4; % Volume of the combustion chamber 

PI_C=9.06e-21*N^4-1.365e-15*N^3+3.3806e-10*N^2-9.7556e-6*N+1.1929;                                                
Cp_a=-1.9747e-028*N^6+ 8.3860e-023*N^5-1.4447e-017*N^4+1.2943e-012*N^3-6.2674e-008*N^2+1.5743e-003*N+9.8846e+002;   
Tt_4=5.5099e-017*N^4-1.1969e-011*N^3+9.7144e-007*N^2-3.7423e-002*N+1.4273e+003;                                      
mf_dot=0.00456; 
A_0=0.001378; 
A_9=0.00302; 

%-----------BETA LINE METHOD---------------------------------------- 
h=0.001;    %Delta "t" 
Tau=0.15;T=0.001;n=1/120000;nce=1/120000;de=0.015;ty=0.5;s=0; 
global nu_p uu 
for j=1:(4/h); 
    N_Ref(j)=85000; 
    e(j)=N_Ref(j)-N(j); 
    if j==1;
        ce(j)=e(j)-0; 
    else  
        ce(j)=e(j)-e(j-1); 
    end 
    en(j)=e(j)*n; 
    cen(j)=ce(j)*nce; 
%-----------------------FUZZYYYY--------- 
%-------RULE 1 -------------------------- 
%if en(i) is P and cen(i) is P then delU_n is P 
    positive(en(j)); 
    nue1=ans; 
    positive(cen(j)); 
    nuce1=ans; 
        if  nue1<=nuce1 
            nu1=nue1; 
        else 
            nu1=nuce1; 
        end 
%-------RULE 2 ---------------------------- 
%if en(i) is N and cen(i) is N then delU_n is N 
    negative(en(j)); 
    nue2=ans; 
    negative(cen(j)); 
    nuce2=ans; 
        if  nue2<=nuce2 
            nu2=nue2; 
        else 
            nu2=nuce2; 
        end 
%-------RULE 3 ---------------------------- 
%if en(i) is P and cen(i) is Z then delU_n is P 
    positive(en(j)); 
    nue3=ans; 
    zero(cen(j)); 
    nuce3=ans; 
        if  nue3<=nuce3 
            nu3=nue3; 
        else 
            nu3=nuce3; 
        end 
%-------RULE 4 ---------------------------- 
%if en(i) is N and cen(i) is Z then delU_n is N 
    negative(en(j)); 
    nue4=ans; 
    zero(cen(j)); 
    nuce4=ans; 
        if  nue4<=nuce4 
            nu4=nue4; 
        else 
            nu4=nuce4; 
        end 
%-------RULE 5 ---------------------------- 
%if en(i) is P and cen(i) is N then delU_n is Z 
    positive(en(j)); 
    nue5=ans; 
    negative(cen(j)); 
    nuce5=ans; 
        if  nue5<=nuce5 
            nu5=nue5; 
        else 
            nu5=nuce5; 
        end 
%-------RULE 6 ---------------------------- 
%if en(i) is N and cen(i) is P then delU_n is Z 
    negative(en(j)); 
    nue6=ans; 
    positive(cen(j)); 
    nuce6=ans; 
    if  nue6<=nuce6 
        nu6=nue6; 
    else 
    nu6=nuce6; 
    end 
%-------RULE 7 ---------------------------- 
%if en(i) is Z and cen(i) is Z then delU_n is Z 
    zero(en(j)); 
    nue7=ans; 
    zero(cen(j)); 
    nuce7=ans; 
    if  nue7<=nuce7 
        nu7=nue7; 
    else 
        nu7=nuce7; 
    end 
%-------RULE 8 ---------------------------- 
%if en(i) is Z and cen(i) is N then delU_n is N 
    zero(en(j)); 
    nue8=ans; 
    negative(cen(j)); 
    nuce8=ans; 
    if  nue8<=nuce8 
        nu8=nue8; 
    else 
        nu8=nuce8; 
    end 
%-------RULE 9 ---------------------------- 
%if en(i) is Z and cen(i) is P then delU_n is P 
    zero(en(j)); 
    nue9=ans; 
    positive(cen(j)); 
    nuce9=ans; 
    if  nue9<=nuce9 
        nu9=nue9; 
    else 
        nu9=nuce9; 
    end 
%-------------------------DEFUZZIFICATION-------------------------- 
dR1=positive2(nu1); 
dR2=negative2(nu2); 
dR3=positive2(nu3); 
dR4=negative2(nu4); 
dR5=zero2(nu5); 
dR6=zero2(nu6); 
dR7=zero2(nu7); 
dR8=negative2(nu8); 
dR9=positive2(nu9); 
d=(nu1*dR1+nu2*dR2+nu3*dR3+nu4*dR4+nu5*dR5+nu6*dR6+nu7*dR7+nu8*dR8+nu9*dR9)/(nu1+nu2+nu3+nu4+nu5+nu6+nu7+nu8+nu9); 
delta_U(j)=d*de; 
if j==1 mf_dot(j)=d*de; 
else     
mf_dot(j)=mf_dot(j-1)+T*d*de; 
end 
%-------------------------BETALINEE----------------------------- 
N_R(j)=N(j)/(teta)^0.5;    
Tt_2(j+1)=Tt_2(j); 
Pt_2(j+1)=Pt_2(j); 
% Corrected RPM 
PI_C_m=[0 1.560 1.648 1.678 0 0 0 0 0 0 0 0  
0 2.019 2.172 2.282 2.337 0 0 0 0 0 0 0   
0 0 2.545 2.750  2.938 3.067 3.121 0 0 0 0 0   
0  0 0 0 3.219 3.438 3.676 3.843 3.955 0 0 0  
0 0 0 0 0 0 3.861 4.090 4.317 4.542 4.695 0]; 
N_cor=[61000 81000 97000 111000 123000]; 
while m==1; 
    if N_cor(i) <= N_R(j) 
        if N_R(j) <= N_cor(i+1) 
            k=(N_cor(i+1)-N_R(j))./(N_cor(i+1)-N_cor(i)); 
            m=0;t=i; 
        end 
    end 
    i=i+1; 
end 
i=1; 
%-------------**----------------% 
for q = 1 : 12; 
N_R_matrix(q)=PI_C_m(t+1,q)-(PI_C_m(t+1,q)-PI_C_m(t,q)).*k; 
end 

while nn==1; 
    if N_R_matrix(c) <= PI_C(j) 
        if PI_C(j) <= N_R_matrix(c+1) 
            Beta=c-(-1*(N_R_matrix(c)-PI_C(j)))./(N_R_matrix(c)-N_R_matrix(c+1)); 
            nn=0; 
        end 
    end 
    c=c+1; 
end 
c=1; 

%-----------Mass Flow Rate Calculations----------------------------- 
x=[1 2 3 4 5 6 7 8 9 10 11 12]; 
y=[61000 81000 97000 111000 123000 ]; 
z=[0 0.2777 0.2185 0.1387 0 0 0 0 0 0 0 0  
0 0.4155 0.3758 0.3226 0.2531 0 0 0 0 0 0 0  
0 0 0.4904 0.4638 0.4358 0.3894 0.318 0 0 0 0 0 
0 0 0 0 0.5206 0.5045 0.4881 0.4514 0.3988 0 0 0  
0 0 0 0 0 0 0.5435 0.5271 0.5094 0.4918 0.4517 0]; 
%-----------Efficiency Calculations--------------------------------- 
x_1=[1 2 3 4 5 6 7 8 9 10 11 12]; 
y_1=[61000 81000 97000 111000 123000 ]; 
z_1=[0 0.711 0.759 0.69 0 0 0 0 0 0 0 0  
0 0.669 0.733 0.76 0.742 0 0 0 0 0 0 0  
0 0 0.646 0.698 0.735 0.752 0.726 0 0 0 0 0  
0 0 0 0 0.65 0.68 0.704 0.715 0.706 0 0 0  
0 0 0 0 0 0 0.62 0.64 0.655 0.664 0.66 0]; 
%------------Results------------------------------------------------ 
m_dot_MAP(j)=interp2(x,y,z,Beta,N_R(j)); 
eff_c(j)=interp2(x_1,y_1,z_1,Beta,N_R(j)); 
%------------------Beta line end------------------------------------ 
%------------------------Compressor--------------------------------- 
mc_dot(j)=m_dot_MAP(j)*(delta)/((teta)^0.5); 
Tt_3(j)=Tt_2(j)*(1+((1/eff_c(j))*((PI_C(j)^((gama_c-1)/gama_c))-1))); 
Pt_3(j)=PI_C(j)*Pt_2(j); 
Ta(j)=(Tt_2(j)+Tt_3(j))/2; 
Cp_a(j+1)=1.0189e3-0.13784*Ta(j)+1.9843e-4*Ta(j)^2+4.2399e-7*Ta(j)^3-3.7632e-10*Ta(j)^4;   
%------------------------Combustion Chamber------------------------- 
Tg(j)=(Tt_4(j)+Tt_3(j))/2; 
Tm(j)=Tg(j); 
f(j)=mf_dot(j)/mc_dot(j); 
Bt(j)=-3.59494e2+4.5164*Tg(j)+2.8116e-3*Tg(j)^2-2.1709e-5*Tg(j)^3+2.8689e-8*Tg(j)^4-1.2263e-11*Tg(j)^5; 
Cp_g(j)=Cp_a(j)+(f(j)/(1+f(j)))*Bt(j); 
%Cp_g(j)=1088.2; 
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
Pt_4_D(j)=Tt_4(j)*R*M_comb(j)/V_comb;  
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
PI_t(j)=-4.5415e-029*N(j)^6+1.7103e-023*N(j)^5-2.6255e-018*N(j)^4+2.1127e-013*N(j)^3-9.4373e-009*N(j)^2+2.1571e-004*N(j)-1.0868e+000; 
Tao_c(j)=Tt_3(j)/Tt_2(j); 
Teta_2(j) =1; 
Teta_3(j) =  Tt_4(j) / Tt_2(j); 
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
%-----------------------THRUST----------- 
Rho_0(j)=-1.4267e-029*N(j)^6+5.8500e-024*N(j)^5-9.7835e-019*N(j)^4+8.5570e-014*N(j)^3-4.1601e-009*N(j)^2+1.0194e-004*N(j)+1.9642e-001; 
m9_dot(j)=mc_dot(j)+ mf_dot(j); 
T_9(j)=8.6554e-026*N(j)^6-3.4961e-020*N(j)^5+5.8197e-015*N(j)^4-5.0939e-010*N(j)^3+2.4663e-005*N(j)^2-6.2630e-001*N(j)+7.2095e+003; 
V9(j)=(gama_t*R*T_9(j))^0.5; 
V0(j)=mc_dot(j)/(Rho_0(j)*A_0); 
F(j)=m9_dot(j)*V9(j)-mc_dot(j)*V0(j); 
end 
subplot(1,3,1);plot(N);grid on;xlabel('t (ms)');ylabel('N (rpm)') 
subplot(1,3,2);plot(Tt_4);grid on;xlabel('t (ms)');ylabel('Tt_4(K)') 
subplot(1,3,3);plot(F);grid on;xlabel('t (ms)');ylabel('Thrust(N)') 


clc; clear; 
format short e 

% Initial parameters and variables
i = 1; m = 1; nn = 1; c = 1; q = 1; v = 1; kl = 0; 
Pt_2 = 93380; % Pascal 
Tt_2(1) = 303;  % Kelvin 
Pt_ref = 100000; % Pascal 
Tt_ref = 298; % Kelvin 
sigma_comb = 0.95;  % Combustion chamber pressure ratio 
nu_comb = 0.6;  % Combustion chamber efficiency
Pi_n = 0.95;  % Nozzle pressure ratio 
M9 = 1;  % Choke condition 
M0 = 0;  % Assumption
H = 0;  % Height
teta = (Tt_2(1) / Tt_ref);  % Dimensionless temperature parameter 
delta = (Pt_2 / Pt_ref);  % Dimensionless pressure parameter 
N(1) = 65000;  % RPM 
I = 1.787e-4;  % Moment of Inertia (kg*m^2) 
eff_t = 0.85;  % Efficiency of the turbine 
gama_c = 1.4;  % Specific heat ratio of compressor 
gama_t = 1.35;  % Specific heat ratio of turbine 
Qf = 43250e3;  % Lower calorific value of fuel [J/kg] 
R = 286.9;  % Individual gas constant [J/(kg*K)] 
V_comb = 6.687233e-4;  % Volume of the combustion chamber

% Compressor and combustion chamber data 
PI_C = 9.06e-21 * N(1)^4 - 1.365e-15 * N(1)^3 + 3.3806e-10 * N(1)^2 - 9.7556e-6 * N(1) + 1.1929;
Cp_a(1) = -1.9747e-028 * N(1)^6 + 8.3860e-023 * N(1)^5 - 1.4447e-017 * N(1)^4 + 1.2943e-012 * N(1)^3 - 6.2674e-008 * N(1)^2 + 1.5743e-003 * N(1) + 9.8846e+002;
Tt_4(1) = 5.5099e-017 * N(1)^4 - 1.1969e-011 * N(1)^3 + 9.7144e-007 * N(1)^2 - 3.7423e-002 * N(1) + 1.4273e+003;
mf_dot(1) = 0.00456;
A_0 = 0.001378;
A_9 = 0.00302;

% Time step
h = 0.001;

% Fuzzy logic control parameters
Tau = 0.15; T = 0.001; n = 1/120000; nce = 1/120000; de = 0.015;
global nu_p uu

% Beta Line Method and Fuzzy Logic Controller
for j = 1:(4/h)
    N_Ref(j) = 85000; 
    e(j) = N_Ref(j) - N(j); 
    
    % Calculate error changes
    if j == 1
        ce(j) = e(j) - 0; 
    else  
        ce(j) = e(j) - e(j-1); 
    end 
    
    en(j) = e(j) * n; 
    cen(j) = ce(j) * nce; 
    
    %---------- Fuzzy Logic Rule Processing --------------------------
    % Combine fuzzy rules using minimums and membership functions
    
    % Rule 1
    positive(en(j)); 
    nue1 = ans; 
    positive(cen(j)); 
    nuce1 = ans; 
    nu1 = min(nue1, nuce1); 
    
    % Rule 2
    negative(en(j)); 
    nue2 = ans; 
    negative(cen(j)); 
    nuce2 = ans; 
    nu2 = min(nue2, nuce2); 
    
    % Rule 3
    positive(en(j)); 
    nue3 = ans; 
    zero(cen(j)); 
    nuce3 = ans; 
    nu3 = min(nue3, nuce3); 
    
    % Rule 4
    negative(en(j)); 
    nue4 = ans; 
    zero(cen(j)); 
    nuce4 = ans; 
    nu4 = min(nue4, nuce4); 
    
    % Rule 5
    positive(en(j)); 
    nue5 = ans; 
    negative(cen(j)); 
    nuce5 = ans; 
    nu5 = min(nue5, nuce5); 
    
    % Rule 6
    negative(en(j)); 
    nue6 = ans; 
    positive(cen(j)); 
    nuce6 = ans; 
    nu6 = min(nue6, nuce6); 
    
    % Rule 7
    zero(en(j)); 
    nue7 = ans; 
    zero(cen(j)); 
    nuce7 = ans; 
    nu7 = min(nue7, nuce7); 
    
    % Rule 8
    zero(en(j)); 
    nue8 = ans; 
    negative(cen(j)); 
    nuce8 = ans; 
    nu8 = min(nue8, nuce8); 
    
    % Rule 9
    zero(en(j)); 
    nue9 = ans; 
    positive(cen(j)); 
    nuce9 = ans; 
    nu9 = min(nue9, nuce9); 
    
    %----------------- Defuzzification -------------------------
    dR1 = positive2(nu1); 
    dR2 = negative2(nu2); 
    dR3 = positive2(nu3); 
    dR4 = negative2(nu4); 
    dR5 = zero2(nu5); 
    dR6 = zero2(nu6); 
    dR7 = zero2(nu7); 
    dR8 = negative2(nu8); 
    dR9 = positive2(nu9); 
    
    % Avoid division by zero during rule aggregation
    total_nu = nu1 + nu2 + nu3 + nu4 + nu5 + nu6 + nu7 + nu8 + nu9;
    if total_nu == 0
        d = 0;
    else
        d = (nu1*dR1 + nu2*dR2 + nu3*dR3 + nu4*dR4 + nu5*dR5 + nu6*dR6 + nu7*dR7 + nu8*dR8 + nu9*dR9) / total_nu;
    end
    
    % Update delta_U
    delta_U(j) = d * de;
    
    if j == 1 
        mf_dot(j) = d * de; 
    else
        mf_dot(j) = mf_dot(j-1) + T * d * de;
    end
    
    %-------------------- BETA LINE METHOD ----------------------------
    N_R(j) = N(j) / (teta)^0.5;
    
    % Corrected RPM compressor map data
    PI_C_m = [
        0 1.560 1.648 1.678 0 0 0 0 0 0 0 0;
        0 2.019 2.172 2.282 2.337 0 0 0 0 0 0 0;
        0 0 2.545 2.750 2.938 3.067 3.121 0 0 0 0 0;
        0 0 0 0 3.219 3.438 3.676 3.843 3.955 0 0 0;
        0 0 0 0 0 0 3.861 4.090 4.317 4.542 4.695 0
    ];
    
    N_cor = [61000 81000 97000 111000 123000];
    
    while m == 1
        if N_cor(i) <= N_R(j)
            if N_R(j) <= N_cor(i+1)
                k = (N_cor(i+1) - N_R(j)) / (N_cor(i+1) - N_cor(i));
                m = 0; t = i;
            end
        end
        i = i + 1;
    end
    i = 1;
    
    % Interpolate between RPM and compressor map data
    for q = 1:12
        N_R_matrix(q) = PI_C_m(t+1, q) - (PI_C_m(t+1, q) - PI_C_m(t, q)) * k;
    end
    
    while nn == 1
        if N_R_matrix(c) <= PI_C(j)
            if PI_C(j) <= N_R_matrix(c+1)
                Beta = c - (-1 * (N_R_matrix(c) - PI_C(j))) / (N_R_matrix(c) - N_R_matrix(c+1));
                nn = 0;
            end
        end
        c = c + 1;
    end
    c = 1;
    
    %-----------Mass Flow Rate and Efficiency Calculations------------
    x = [1 2 3 4 5 6 7 8 9 10 11 12];
    y = [61000 81000 97000 111000 123000];
    z = [
        0 0.2777 0.2185 0.1387 0 0 0 0 0 0 0 0;
        0 0.4155 0.3758 0.3226 0.2531 0 0 0 0 0 0 0;
        0 0 0.4904 0.4638 0.4358 0.3894 0.318 0 0 0 0 0;
        0 0 0 0 0.5206 0.5045 0.4881 0.4514 0.3988 0 0 0;
        0 0 0 0 0 0 0.5435 0.5271 0.5094 0.4918 0.4517 0
    ];
    
    x_1 = [1 2 3 4 5 6 7 8 9 10 11 12];
    y_1 = [61000 81000 97000 111000 123000];
    z_1 = [
        0 0.711 0.759 0.69 0 0 0 0 0 0 0 0;
        0 0.669 0.733 0.76 0.742 0 0 0 0 0 0 0;
        0 0 0.646 0.698 0.735 0.752 0.726 0 0 0 0 0;
        0 0 0 0 0.65 0.68 0.704 0.715 0.706 0 0 0;
        0 0 0 0 0 0 0.62 0.64 0.655 0.664 0.66 0
    ];
    
    % Get mass flow rate and efficiency from interpolation
    m_dot_MAP(j) = interp2(x, y, z, Beta, N_R(j));
    eff_c(j) = interp2(x_1, y_1, z_1, Beta, N_R(j));
    
    %------------------------Compressor--------------------------------- 
    mc_dot(j) = m_dot_MAP(j) * (delta) / ((teta)^0.5); 
    Tt_3(j) = Tt_2(j) * (1 + ((1 / eff_c(j)) * ((PI_C(j)^((gama_c-1) / gama_c)) - 1)));
    Pt_3(j) = PI_C(j) * Pt_2(j);
    Ta(j) = (Tt_2(j) + Tt_3(j)) / 2;
    Cp_a(j+1) = 1.0189e3 - 0.13784 * Ta(j) + 1.9843e-4 * Ta(j)^2 + 4.2399e-7 * Ta(j)^3 - 3.7632e-10 * Ta(j)^4;
    
    %------------------------Combustion Chamber-------------------------
    Tg(j) = (Tt_4(j) + Tt_3(j)) / 2;
    Tm(j) = Tg(j);
    f(j) = mf_dot(j) / mc_dot(j);
    Bt(j) = -3.59494e2 + 4.5164 * Tg(j) + 2.8116e-3 * Tg(j)^2 - 2.1709e-5 * Tg(j)^3 + 2.8689e-8 * Tg(j)^4 - 1.2263e-11 * Tg(j)^5;
    Cp_g(j) = Cp_a(j) + (f(j) / (1 + f(j))) * Bt(j);
    Cv_comb(j) = Cp_g(j) - R;
    
    Pt_4(j) = sigma_comb * Pt_3(j);
    Pm(j) = (Pt_4(j) + Pt_3(j)) / 2;
    Rho(j) = Pm(j) / (R * Tm(j));
    M_comb(j) = Rho(j) * V_comb;
    
    %------------------Runge Kutta fourth order for Tt_4 ---------------  
    F(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j)); 
    K1(j) = F(j);
    K2(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + (h/2) * K1(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j)); 
    K3(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + (h/2) * K2(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j));
    K4(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + h * K3(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j));
    Tt_4(j+1) = Tt_4(j) + (h/6) * (K1(j) + 2 * K2(j) + 2 * K3(j) + K4(j));
    
    %------------------Runge Kutta fourth order for Pt_4----------------  
    H(j) = Pt_4(j) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j)); 
    M1(j) = H(j); 
    M2(j) = (Pt_4(j) + h/2 * M1(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    M3(j) = (Pt_4(j) + h/2 * M2(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    M4(j) = (Pt_4(j) + h * M3(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    Pt_4(j+1) = Pt_4(j) + h/6 * (M1(j) + 2 * M2(j) + 2 * M3(j) + M4(j));
    
    Pt_3(j+1) = Pt_4(j+1) / sigma_comb; 
    PI_C(j+1) = Pt_3(j+1) / Pt_2(j);
    
    %--------------------TURBINE---------------------------------------- 
    PI_t(j) = -4.5415e-029 * N(j)^6 + 1.7103e-023 * N(j)^5 - 2.6255e-018 * N(j)^4 + 2.1127e-013 * N(j)^3 - 9.4373e-009 * N(j)^2 + 2.1571e-004 * N(j) - 1.0868;
    Tao_c(j) = Tt_3(j) / Tt_2(j); 
    Teta_2(j) = 1;
    Teta_3(j) = Tt_4(j) / Tt_2(j);
    Tao_t(j) = 1 - ((Cp_a(j) * Teta_2(j)) / (Cp_g(j) * (1 + f(j)) * Teta_3(j))) * (Tao_c(j) - 1);
    Tt_5(j) = 685;
    
    %--------------------Work Balance----------------------------------- 
    Pt(j) = ((mc_dot(j) + mf_dot(j)) * (Tt_4(j) - Tt_5(j)) * Cp_g(j));
    Pc(j) = mc_dot(j) * Cp_a(j) * (Tt_3(j) - Tt_2(j)); 
    Pf(j) = 500; 
    if Pf(j) < 0
        Pf(j) = 0;
    end
    W(j) = N(j) * pi / 30; 
    W(j+1) = sqrt((2 * (Pt(j) - Pc(j) - Pf(j)) / I) * h + W(j)^2); 
    N(j+1) = W(j+1) * 30 / pi;
    
    %-----------------------THRUST------------------------------------- 
    Rho_0(j) = -1.4267e-029 * N(j)^6 + 5.8500e-024 * N(j)^5 - 9.7835e-019 * N(j)^4 + 8.5570e-014 * N(j)^3 - 4.1601e-009 * N(j)^2 + 1.0194e-004 * N(j) + 1.9642e-001; 
    m9_dot(j) = mc_dot(j) + mf_dot(j); 
    T_9(j) = 8.6554e-026 * N(j)^6 - 3.4961e-020 * N(j)^5 + 5.8197e-015 * N(j)^4 - 5.0939e-010 * N(j)^3 + 2.4663e-005 * N(j)^2 - 6.2630e-001 * N(j) + 7.2095e+003;
    V9(j) = (gama_t * R * T_9(j))^0.5; 
    V0(j) = mc_dot(j) / (Rho_0(j) * A_0); 
    F(j) = m9_dot(j) * V9(j) - mc_dot(j) * V0(j);
end

%-----------------------PLOT RESULTS---------------------------------
subplot(1,3,1); plot(N); grid on; xlabel('t (ms)'); ylabel('N (rpm)'); 
subplot(1,3,2); plot(Tt_4); grid on; xlabel('t (ms)'); ylabel('Tt_4 (K)'); 
subplot(1,3,3); plot(F); grid on; xlabel('t (ms)'); ylabel('Thrust (N)');

%}

clc; clear; 
format short e 

% Initial parameters and variables
i = 1; m = 1; nn = 1; c = 1; q = 1; v = 1; kl = 0; 
Pt_2(1) = 93380; % Pascal (initialize as array)
Tt_2(1) = 303;  % Kelvin (initialize as array)
Pt_ref = 100000; % Pascal 
Tt_ref = 298; % Kelvin 
sigma_comb = 0.95;  % Combustion chamber pressure ratio 
nu_comb = 0.6;  % Combustion chamber efficiency
Pi_n = 0.95;  % Nozzle pressure ratio 
M9 = 1;  % Choke condition 
M0 = 0;  % Assumption
H = 0;  % Height
teta = (Tt_2(1) / Tt_ref);  % Dimensionless temperature parameter 
delta = (Pt_2(1) / Pt_ref);  % Dimensionless pressure parameter 
N(1) = 65000;  % RPM (initialize as array)
I = 1.787e-4;  % Moment of Inertia (kg*m^2) 
eff_t = 0.85;  % Efficiency of the turbine 
gama_c = 1.4;  % Specific heat ratio of compressor 
gama_t = 1.35;  % Specific heat ratio of turbine 
Qf = 43250e3;  % Lower calorific value of fuel [J/kg] 
R = 286.9;  % Individual gas constant [J/(kg*K)] 
V_comb = 6.687233e-4;  % Volume of the combustion chamber

% Compressor and combustion chamber data 
PI_C(1) = 9.06e-21 * N(1)^4 - 1.365e-15 * N(1)^3 + 3.3806e-10 * N(1)^2 - 9.7556e-6 * N(1) + 1.1929;
Cp_a(1) = -1.9747e-028 * N(1)^6 + 8.3860e-023 * N(1)^5 - 1.4447e-017 * N(1)^4 + 1.2943e-012 * N(1)^3 - 6.2674e-008 * N(1)^2 + 1.5743e-003 * N(1) + 9.8846e+002;
Tt_4(1) = 5.5099e-017 * N(1)^4 - 1.1969e-011 * N(1)^3 + 9.7144e-007 * N(1)^2 - 3.7423e-002 * N(1) + 1.4273e+003;
mf_dot(1) = 0.00456;
A_0 = 0.001378;
A_9 = 0.00302;

% Time step
h = 0.001;

% Fuzzy logic control parameters
Tau = 0.15; T = 0.001; n = 1/120000; nce = 1/120000; de = 0.015;
global nu_p uu

% Initialize m_dot_MAP, eff_c, Tt_2, Pt_2 as arrays to hold values
m_dot_MAP = zeros(1, (4/h)); 
eff_c = zeros(1, (4/h));
Tt_2 = zeros(1, (4/h)); Tt_2(1) = 303; % initialize first value
Pt_2 = zeros(1, (4/h)); Pt_2(1) = 93380; % initialize first value

% Beta Line Method and Fuzzy Logic Controller
for j = 1:(4/h)
    N_Ref(j) = 85000; 
    e(j) = N_Ref(j) - N(j); 
    
    % Calculate error changes
    if j == 1
        ce(j) = e(j) - 0; 
    else  
        ce(j) = e(j) - e(j-1); 
    end 
    
    en(j) = e(j) * n; 
    cen(j) = ce(j) * nce; 
    
    %---------- Fuzzy Logic Rule Processing --------------------------
    % Combine fuzzy rules using minimums and membership functions
    % Rule 1
    positive(en(j)); 
    nue1 = ans; 
    positive(cen(j)); 
    nuce1 = ans; 
    nu1 = min(nue1, nuce1); 
    
    % Rule 2
    negative(en(j)); 
    nue2 = ans; 
    negative(cen(j)); 
    nuce2 = ans; 
    nu2 = min(nue2, nuce2); 
    
    % Rule 3
    positive(en(j)); 
    nue3 = ans; 
    zero(cen(j)); 
    nuce3 = ans; 
    nu3 = min(nue3, nuce3); 
    
    % Rule 4
    negative(en(j)); 
    nue4 = ans; 
    zero(cen(j)); 
    nuce4 = ans; 
    nu4 = min(nue4, nuce4); 
    
    % Rule 5
    positive(en(j)); 
    nue5 = ans; 
    negative(cen(j)); 
    nuce5 = ans; 
    nu5 = min(nue5, nuce5); 
    
    % Rule 6
    negative(en(j)); 
    nue6 = ans; 
    positive(cen(j)); 
    nuce6 = ans; 
    nu6 = min(nue6, nuce6); 
    
    % Rule 7
    zero(en(j)); 
    nue7 = ans; 
    zero(cen(j)); 
    nuce7 = ans; 
    nu7 = min(nue7, nuce7); 
    
    % Rule 8
    zero(en(j)); 
    nue8 = ans; 
    negative(cen(j)); 
    nuce8 = ans; 
    nu8 = min(nue8, nuce8); 
    
    % Rule 9
    zero(en(j)); 
    nue9 = ans; 
    positive(cen(j)); 
    nuce9 = ans; 
    nu9 = min(nue9, nuce9); 
    
    %----------------- Defuzzification -------------------------
    dR1 = positive2(nu1); 
    dR2 = negative2(nu2); 
    dR3 = positive2(nu3); 
    dR4 = negative2(nu4); 
    dR5 = zero2(nu5); 
    dR6 = zero2(nu6); 
    dR7 = zero2(nu7); 
    dR8 = negative2(nu8); 
    dR9 = positive2(nu9); 
    
    % Avoid division by zero during rule aggregation
    total_nu = nu1 + nu2 + nu3 + nu4 + nu5 + nu6 + nu7 + nu8 + nu9;
    if total_nu == 0
        d = 0;
    else
        d = (nu1*dR1 + nu2*dR2 + nu3*dR3 + nu4*dR4 + nu5*dR5 + nu6*dR6 + nu7*dR7 + nu8*dR8 + nu9*dR9) / total_nu;
    end
    
    % Update delta_U
    delta_U(j) = d * de;
    
    if j == 1 
        mf_dot(j) = d * de; 
    else
        mf_dot(j) = mf_dot(j-1) + T * d * de;
    end
    
    %-------------------- BETA LINE METHOD ----------------------------
    N_R(j) = N(j) / (teta)^0.5;
    Tt_2(j+1)=Tt_2(j); 
    Pt_2(j+1)=Pt_2(j); 
    
    % Corrected RPM compressor map data
    PI_C_m = [
        0 1.560 1.648 1.678 0 0 0 0 0 0 0 0;
        0 2.019 2.172 2.282 2.337 0 0 0 0 0 0 0;
        0 0 2.545 2.750 2.938 3.067 3.121 0 0 0 0 0;
        0 0 0 0 3.219 3.438 3.676 3.843 3.955 0 0 0;
         0 0 0 0 0 0 3.861 4.090 4.317 4.542 4.695 0
    ];
    
    N_cor = [61000 81000 97000 111000 123000];
    
    while m == 1
        if N_cor(i) <= N_R(j)
            if N_R(j) <= N_cor(i+1)
                k = (N_cor(i+1) - N_R(j)) / (N_cor(i+1) - N_cor(i));
                m = 0; t = i;
            end
        end
        i = i + 1;
    end
    i = 1;
    
    % Interpolate between RPM and compressor map data
    for q = 1:12
        N_R_matrix(q) = PI_C_m(t+1, q) - (PI_C_m(t+1, q) - PI_C_m(t, q)) * k;
    end
    
    while nn == 1
        if N_R_matrix(c) <= PI_C(j)
            if PI_C(j) <= N_R_matrix(c+1)
                Beta = c - (-1 * (N_R_matrix(c) - PI_C(j))) / (N_R_matrix(c) - N_R_matrix(c+1));
                nn = 0;
            end
        end
        c = c + 1;
    end
    c = 1;
    
    %-----------Mass Flow Rate and Efficiency Calculations------------
    x = [1 2 3 4 5 6 7 8 9 10 11 12];
    y = [61000 81000 97000 111000 123000];
    z = [
        0 0.2777 0.2185 0.1387 0 0 0 0 0 0 0 0;
        0 0.4155 0.3758 0.3226 0.2531 0 0 0 0 0 0 0;
        0 0 0.4904 0.4638 0.4358 0.3894 0.318 0 0 0 0 0;
        0 0 0 0 0.5206 0.5045 0.4881 0.4514 0.3988 0 0 0;
        0 0 0 0 0 0 0.5435 0.5271 0.5094 0.4918 0.4517 0
    ];
    
    x_1 = [1 2 3 4 5 6 7 8 9 10 11 12];
    y_1 = [61000 81000 97000 111000 123000];
    z_1 = [
        0 0.711 0.759 0.69 0 0 0 0 0 0 0 0;
        0 0.669 0.733 0.76 0.742 0 0 0 0 0 0 0;
        0 0 0.646 0.698 0.735 0.752 0.726 0 0 0 0 0;
        0 0 0 0 0.65 0.68 0.704 0.715 0.706 0 0 0;
        0 0 0 0 0 0 0.62 0.64 0.655 0.664 0.66 0
    ];
    
    % Get mass flow rate and efficiency from interpolation
    m_dot_MAP(j) = interp2(x, y, z, Beta, N_R(j));  % Interpolated mass flow rate
    eff_c(j) = interp2(x_1, y_1, z_1, Beta, N_R(j));  % Interpolated efficiency
    
    %------------------------Compressor--------------------------------- 
    mc_dot(j) = m_dot_MAP(j) * (delta) / ((teta)^0.5); 
    Tt_3(j) = Tt_2(j) * (1 + ((1 / eff_c(j)) * ((PI_C(j)^((gama_c-1) / gama_c)) - 1)));
    Pt_3(j) = PI_C(j) * Pt_2(j);
    Ta(j) = (Tt_2(j) + Tt_3(j)) / 2;
    Cp_a(j+1) = 1.0189e3 - 0.13784 * Ta(j) + 1.9843e-4 * Ta(j)^2 + 4.2399e-7 * Ta(j)^3 - 3.7632e-10 * Ta(j)^4;
    
    %------------------------Combustion Chamber-------------------------
    Tg(j) = (Tt_4(j) + Tt_3(j)) / 2;
    Tm(j) = Tg(j);
    f(j) = mf_dot(j) / mc_dot(j);
    Bt(j) = -3.59494e2 + 4.5164 * Tg(j) + 2.8116e-3 * Tg(j)^2 - 2.1709e-5 * Tg(j)^3 + 2.8689e-8 * Tg(j)^4 - 1.2263e-11 * Tg(j)^5;
    Cp_g(j) = Cp_a(j) + (f(j) / (1 + f(j))) * Bt(j);
    Cv_comb(j) = Cp_g(j) - R;
    
    Pt_4(j) = sigma_comb * Pt_3(j);
    Pm(j) = (Pt_4(j) + Pt_3(j)) / 2;
    Rho(j) = Pm(j) / (R * Tm(j));
    M_comb(j) = Rho(j) * V_comb;
    
    %------------------Runge Kutta fourth order for Tt_4 ---------------  
    F(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j)); 
    K1(j) = F(j);
    K2(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + (h/2) * K1(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j)); 
    K3(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + (h/2) * K2(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j));
    K4(j) = (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * (Tt_4(j) + h * K3(j)) + Qf * nu_comb * mf_dot(j)) / (Cv_comb(j) * M_comb(j));
    Tt_4(j+1) = Tt_4(j) + (h/6) * (K1(j) + 2 * K2(j) + 2 * K3(j) + K4(j));
    
    %------------------Runge Kutta fourth order for Pt_4----------------  
    H(j) = Pt_4(j) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j)); 
    M1(j) = H(j); 
    M2(j) = (Pt_4(j) + h/2 * M1(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    M3(j) = (Pt_4(j) + h/2 * M2(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    M4(j) = (Pt_4(j) + h * M3(j)) * (mc_dot(j) * Cp_a(j) * Tt_3(j) - (mc_dot(j) + mf_dot(j)) * Cp_g(j) * Tt_4(j) + Qf * nu_comb * mf_dot(j)) / (Tt_4(j) * Cv_comb(j) * M_comb(j));
    Pt_4(j+1) = Pt_4(j) + h/6 * (M1(j) + 2 * M2(j) + 2 * M3(j) + M4(j));
    
    Pt_3(j+1) = Pt_4(j+1) / sigma_comb; 
    PI_C(j+1) = Pt_3(j+1) / Pt_2(j);
    
    %--------------------TURBINE---------------------------------------- 
    PI_t(j) = -4.5415e-029 * N(j)^6 + 1.7103e-023 * N(j)^5 - 2.6255e-018 * N(j)^4 + 2.1127e-013 * N(j)^3 - 9.4373e-009 * N(j)^2 + 2.1571e-004 * N(j) - 1.0868;
    Tao_c(j) = Tt_3(j) / Tt_2(j); 
    Teta_2(j) = 1;
    Teta_3(j) = Tt_4(j) / Tt_2(j);
    Tao_t(j) = 1 - ((Cp_a(j) * Teta_2(j)) / (Cp_g(j) * (1 + f(j)) * Teta_3(j))) * (Tao_c(j) - 1);
    Tt_5(j) = 685;
    
    %--------------------Work Balance----------------------------------- 
    Pt(j) = ((mc_dot(j) + mf_dot(j)) * (Tt_4(j) - Tt_5(j)) * Cp_g(j));
    Pc(j) = mc_dot(j) * Cp_a(j) * (Tt_3(j) - Tt_2(j)); 
    Pf(j) = 500; 
    if Pf(j) < 0
        Pf(j) = 0;
    end
    W(j) = N(j) * pi / 30; 
    W(j+1) = sqrt((2 * (Pt(j) - Pc(j) - Pf(j)) / I) * h + W(j)^2); 
    N(j+1) = W(j+1) * 30 / pi;
    
    %-----------------------THRUST------------------------------------- 
    Rho_0(j) = -1.4267e-029 * N(j)^6 + 5.8500e-024 * N(j)^5 - 9.7835e-019 * N(j)^4 + 8.5570e-014 * N(j)^3 - 4.1601e-009 * N(j)^2 + 1.0194e-004 * N(j) + 1.9642e-001; 
    m9_dot(j) = mc_dot(j) + mf_dot(j); 
    T_9(j) = 8.6554e-026 * N(j)^6 - 3.4961e-020 * N(j)^5 + 5.8197e-015 * N(j)^4 - 5.0939e-010 * N(j)^3 + 2.4663e-005 * N(j)^2 - 6.2630e-001 * N(j) + 7.2095e+003;
    V9(j) = (gama_t * R * T_9(j))^0.5; 
    V0(j) = mc_dot(j) / (Rho_0(j) * A_0); 
    F(j) = m9_dot(j) * V9(j) - mc_dot(j) * V0(j);
end

%-----------------------PLOT RESULTS---------------------------------
subplot(1,3,1); plot(N,'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('N (rpm)'); 
subplot(1,3,2); plot(Tt_4,'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('Tt_4 (K)'); 
subplot(1,3,3); plot(F,'LineWidth', 2); grid on; xlabel('t (ms)'); ylabel('Thrust (N)');

