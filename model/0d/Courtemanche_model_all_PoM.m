function [dY, ionCurrents] = Courtemanche_model_all_PoM(time_instat, Y, ~, settings, param_factors)

%-----------------------------------------------------------------------------
% Definition of state variables
%-----------------------------------------------------------------------------

i_u      = 1;  i_v      = 2;  i_w      = 3;  i_d      = 4;   i_f_Ca   = 5;  i_f      = 6; 
i_h      = 7;  i_j      = 8;  i_m      = 9;  i_Ca_i   = 10;  i_Ca_rel = 11; i_Ca_up  = 12; 
i_K_i    = 13; i_Na_i   = 14; i_V      = 15; i_xr     = 16;  i_xs     = 17; i_oa     = 18; 
i_oi     = 19; i_ua     = 20; i_ui     = 21; i_f_Ca_sk = 22; i_O_K2P = 23; 

dY = zeros(length(Y),1);


%-----------------------------------------------------------------------------
% Select model variant
%-----------------------------------------------------------------------------

atrial_model = settings.atrial_model;
AF_condition = settings.AF_condition;

%%% atrial model (electrical heterogeneity)
switch atrial_model
    case 'RA_PM'
        fact_het_Ito  = 1;          fact_het_ICaL = 1;          fact_het_IKr  = 1;  
    case 'CT_BBra' 
        fact_het_Ito  = 1;          fact_het_ICaL = 1.67;       fact_het_IKr  = 1;  
    case 'TVR' 
        fact_het_Ito  = 1;          fact_het_ICaL = 0.67;       fact_het_IKr  = 1.53;
    case 'RAA' 
        fact_het_Ito  = 0.68;       fact_het_ICaL = 1;          fact_het_IKr  = 1;
    case 'LA' 
        fact_het_Ito  = 1;          fact_het_ICaL = 1;          fact_het_IKr  = 1.6;  
    case 'BBla' 
        fact_het_Ito  = 1;          fact_het_ICaL = 1.67;       fact_het_IKr  = 1.6;
    case 'PV' 
        fact_het_Ito  = 1;          fact_het_ICaL = 1;          fact_het_IKr  = 2.2;
    case 'LAA' 
        fact_het_Ito  = 0.68;       fact_het_ICaL = 1;          fact_het_IKr  = 1.6;
    case 'MVR' 
        fact_het_Ito  = 1;          fact_het_ICaL = 0.67;       fact_het_IKr  = 2.44;
    otherwise
        disp('\n Unknown atrial model')
end

%%% AF condition used in the manuscript
switch AF_condition
    case 'pAF'
        % Paper pAF population (called nSR in the legacy simulation folders)
        fact_AF_Ito = 1; fact_AF_ICaL = 1; fact_AF_IK1 = 1; fact_AF_IKur = 1;
        fact_AF_IKs = 1; fact_AF_IKACh = 1; fact_AF_IKr = 1; fact_AF_INa = 1;
    case 'IREpAF'
        % Paper IREpAF population (called pAF in the legacy simulation folders)
        fact_AF_Ito = 1; fact_AF_ICaL = 1; fact_AF_IK1 = 2.00; fact_AF_IKur = 1;
        fact_AF_IKs = 1; fact_AF_IKACh = 2; fact_AF_IKr = 1; fact_AF_INa = 1;
    otherwise
        error('Unknown AF condition: %s. Expected pAF or IREpAF.', AF_condition);
end

%%% DRUG application (pore simple model)
D_drug = settings.drug_conc;        % Block = 1/(1 + (D_dof/IC50)^nH)
fact_drug_ICaL  = 1/(1 + ( D_drug / settings.drug_ICaL_IC50  )^settings.drug_ICaL_nH );
fact_drug_IK1   = 1/(1 + ( D_drug / settings.drug_IK1_IC50   )^settings.drug_IK1_nH  );
fact_drug_IKACh = 1/(1 + ( D_drug / settings.drug_IKACh_IC50 )^settings.drug_IKACh_nH);
fact_drug_IKr   = 1/(1 + ( D_drug / settings.drug_IKr_IC50   )^settings.drug_IKr_nH  );
fact_drug_IKs   = 1/(1 + ( D_drug / settings.drug_IKs_IC50   )^settings.drug_IKs_nH  );
fact_drug_IKur  = 1/(1 + ( D_drug / settings.drug_IKur_IC50  )^settings.drug_IKur_nH );
fact_drug_ISK   = 1/(1 + ( D_drug / settings.drug_ISK_IC50   )^settings.drug_ISK_nH  );
fact_drug_IK2P  = 1/(1 + ( D_drug / settings.drug_IK2P_IC50  )^settings.drug_IK2P_nH  );
fact_drug_Ito   = 1/(1 + ( D_drug / settings.drug_Ito_IC50   )^settings.drug_Ito_nH  );
fact_drug_INa   = 1/(1 + ( D_drug / settings.drug_INa_IC50   )^settings.drug_INa_nH  );
fact_drug_INaK  = 1/(1 + ( D_drug / settings.drug_INaK_IC50  )^settings.drug_INaK_nH );
fact_drug_INCX  = 1/(1 + ( D_drug / settings.drug_INCX_IC50  )^settings.drug_INCX_nH );


fact_curr_ICaL  = param_factors(1);
fact_curr_IK1   = param_factors(2);
fact_curr_IKr   = param_factors(3);
fact_curr_IKs   = param_factors(4);
fact_curr_IKur  = param_factors(5);
fact_curr_Ito   = param_factors(6);
fact_curr_IKACh = param_factors(7);
fact_curr_INa   = param_factors(8);
fact_curr_INaK  = param_factors(9);
fact_curr_INCX  = param_factors(10);
fact_curr_ISK   = param_factors(11);
fact_curr_IK2P  = param_factors(12);
fact_curr_ISAC  = param_factors(13);

%%% total factor  =        currF          *      het       *      AF       *      DRUG
fact_tot_ICaL  =          fact_curr_ICaL  * fact_het_ICaL  * fact_AF_ICaL  * fact_drug_ICaL;
fact_tot_IK1   =          fact_curr_IK1                    * fact_AF_IK1   * fact_drug_IK1;
fact_tot_IKACh =          fact_curr_IKACh                  * fact_AF_IKACh * fact_drug_IKACh;
fact_tot_IKr   =          fact_curr_IKr   * fact_het_IKr   * fact_AF_IKr   * fact_drug_IKr;
fact_tot_IKs   =          fact_curr_IKs                    * fact_AF_IKs   * fact_drug_IKs;
fact_tot_IKur  =          fact_curr_IKur                   * fact_AF_IKur  * fact_drug_IKur;
fact_tot_Ito   =          fact_curr_Ito   * fact_het_Ito   * fact_AF_Ito   * fact_drug_Ito;
fact_tot_ISK   =          fact_curr_ISK                                    * fact_drug_ISK;
fact_tot_IK2P  =          fact_curr_IK2P                                   * fact_drug_IK2P;
fact_tot_ISAC  =          fact_curr_ISAC;
fact_tot_INa   =          fact_curr_INa                    * fact_AF_INa   * fact_drug_INa;
fact_tot_INaK  =          fact_curr_INaK                                   * fact_drug_INaK;
fact_tot_INCX  =          fact_curr_INCX                                   * fact_drug_INCX;
fact_tot_Irel  = settings.fact_curr_Irel; 
fact_tot_Ileak = settings.fact_curr_Ileak;
fact_tot_Iup   = settings.fact_curr_Iup; 
fact_tot_Itr   = settings.fact_curr_Itr; 
fact_tot_ICaP  = settings.fact_curr_ICaP; 


%-----------------------------------------------------------------------------
% Model parameters
%-----------------------------------------------------------------------------

% Universal constants
F = 96.4867;                % coulomb_per_millimole (in membrane)
R = 8.3143;                 % joule_per_mole_kelvin (in membrane)
T = 310.0;                  % kelvin (in membrane)

% Cell constants
Cm = 100.0;                 % picoF (in membrane) 
V_cell = 20100.0;           % micrometre_3 (in intracellular_ion_concentrations)
V_i = V_cell*0.68;          % micrometre_3 (in intracellular_ion_concentrations)
V_up = V_cell*0.0552;       % micrometre_3 (in intracellular_ion_concentrations)
V_rel = V_cell*0.0048;      % micrometre_3 (in intracellular_ion_concentrations)

% Channel constants
CMDN_max = 0.05;            % millimolar (in Ca_buffers)
CSQN_max = 10.0;            % millimolar (in Ca_buffers)
Km_CMDN = 0.00238;          % millimolar (in Ca_buffers)
Km_CSQN = 0.8;              % millimolar (in Ca_buffers)
Km_TRPN = 0.0005;           % millimolar (in Ca_buffers)
TRPN_max = 0.07;            % millimolar (in Ca_buffers)
Ca_up_max = 15.0;           % millimolar (in Ca_leak_current_by_the_NSR)
K_rel = 30.0;               % per_millisecond -1 (in Ca_release_current_from_JSR)
I_up_max = 0.005;           % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
K_up = 0.00092;             % millimolar (in Ca_uptake_current_by_the_NSR)
I_NaCa_max = 1600.0;        % picoA_per_picoF (in Na_Ca_exchanger_current)
K_mCa = 1.38;               % millimolar (in Na_Ca_exchanger_current)
K_mNa = 87.5;               % millimolar (in Na_Ca_exchanger_current)
K_sat = 0.1;                % dimensionless (in Na_Ca_exchanger_current)
gamma = 0.35;               % dimensionless (in Na_Ca_exchanger_current)
g_B_Ca = 0.001131;          % nanoS_per_picoF (in background_currents)
g_B_K = 0.0;                % nanoS_per_picoF (in background_currents)
g_B_Na = 0.0006744375;      % nanoS_per_picoF (in background_currents)
g_Kr = .0294;               % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
i_CaP_max = 0.275;          % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
g_Ks = 0.12941176;          % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
Km_K_o = 1.5;               % millimolar (in sodium_potassium_pump)
Km_Na_i = 10.0;             % millimolar (in sodium_potassium_pump)
i_NaK_max = 0.59933874;     % picoA_per_picoF (in sodium_potassium_pump)
Ca_o = 1.8;                 % millimolar (in standard_ionic_concentrations)
K_o = 5.4;                  % millimolar (in standard_ionic_concentrations)
Na_o = 140.0;               % millimolar (in standard_ionic_concentrations)
g_K1 = 0.09;                % nanoS_per_picoF (in time_independent_potassium_current)
tau_tr = 180.0;             % millisecond (in transfer_current_from_NSR_to_JSR) %equation n�70
K_Q10 = 3.0;                % dimensionless (in transient_outward_K_current)
g_to = 0.1652;              % nanoS_per_picoF (in transient_outward_K_current)
g_Ca_L = 0.12375;           % nanoS_per_picoF (in L_type_Ca_channel)
g_Na = 7.8;                 % nanoS_per_picoF (in Na_channel)

% SK and K2P conductivities
g_SK  = settings.gSK;       % conductivity of SK channel [nS/pF]
g_K2P = settings.gK2P;      % conductivity of K2P channel [nS/pF]

% SAC parameters
lambda = settings.SAC_lambda;   % strech ratio (1 if cell is not stretched)
K_sac = 100;                    % Parameter to define Isac when not stretched [-]
alpha_sac = 3;                  % Parameter to describe sensitivity to stretch [-]
P_Na = 1;                       % Channel relative permeability to Na+
P_Ca = 1;                       % Channel relative permeability to Ca2+
P_K  = 1;                       % Channel relative permeability to K+
Z_Na = 1;                       % Na+ ion valence
Z_Ca = 2;                       % Ca2+ ion valence
Z_K  = 1;                       % K+ ion valence
G_sac = settings.gSAC;          % 0.05 - maximun membrane conductance
g_sac = G_sac/(1 + K_sac*exp(-alpha_sac*(lambda-1))); % conductivity of SAC channel [nS]


%-----------------------------------------------------------------------------
% Stimulus current (modificado Violeta)
%-----------------------------------------------------------------------------
if (time_instat < settings.stim_dur+settings.stim_offset)  && (time_instat >= settings.stim_offset) 
    i_st = settings.stim_amp;
else
    i_st = 0.0;
end


%-----------------------------------------------------------------------------
% Ion Currents 
%-----------------------------------------------------------------------------

%%% Equilibrium potentials
E_Na = R*T/F*log(Na_o/Y(i_Na_i)); %(Eq.28)
E_Ca = R*T/(2.0*F)*log(Ca_o/Y(i_Ca_i)); %(Eq.28)
E_K = R*T/F*log(K_o/Y(i_K_i)); %(Eq.28)

%%% Fast Na+ Current       
INa = 1.3*Cm*fact_tot_INa*g_Na*Y(i_m)^3.0*Y(i_h)*Y(i_j)*(Y(i_V)-E_Na); %(Eq.29)

% m
if (Y(i_V) == -47.13) %(Eq.30)
   alpha_m = 3.2;
else
   alpha_m = 0.32*(Y(i_V)+47.13)/(1.0-exp(-0.1*(Y(i_V)+47.13)));
end
beta_m = 0.08*exp(-Y(i_V)/11.0); %(Eq.30)
m_inf = alpha_m/(alpha_m+beta_m); %(Eq.34)
tau_m = (alpha_m+beta_m)^(-1); %(Eq.34)

% h
if (Y(i_V) < -40.0) %(Eq.31)
   alpha_h = 0.135*exp((Y(i_V)+80.0)/-6.8);
else
   alpha_h = 0.0;
end
if (Y(i_V) < -40.0) %(Eq.31)
   beta_h = 3.56*exp(0.079*Y(i_V))+3.1e5*exp(0.35*Y(i_V));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(i_V)+10.66)/-11.1)));
end
h_inf = alpha_h/(alpha_h+beta_h);
tau_h = 1.0/(alpha_h+beta_h); %(Eq.34)

% j 
if (Y(i_V) < -40.0) %(Eq.32)
   alpha_j = (-1.2714e5*exp(0.2444*Y(i_V))-3.474e-5*exp(-0.04391*Y(i_V)))*(Y(i_V)+37.78)/(1.0+exp(0.311*(Y(i_V)+79.23)));
else
   alpha_j = 0.0;
end
if (Y(i_V) < -40.0) %(Eq.33)
   beta_j = 0.1212*exp(-0.01052*Y(i_V))/(1.0+exp(-0.1378*(Y(i_V)+40.14)));
else
   beta_j = 0.3*exp(-2.535e-7*Y(i_V))/(1.0+exp(-0.1*(Y(i_V)+32.0)));
end
j_inf = alpha_j/(alpha_j+beta_j); %(Eq.34)
tau_j = 1.0/(alpha_j+beta_j); %(Eq.34)


%%% Time-Independent K+ Current
IK1 = Cm*fact_tot_IK1*g_K1*(Y(i_V)-E_K)/(1.0+exp(0.07*(Y(i_V)+80.0))); %(Eq.35)


%%% Transient Outward K+ Current
Ito = Cm*fact_tot_Ito*g_to*Y(i_oa)^3.0*Y(i_oi)*(Y(i_V)-E_K); %(Eq.36)

% oa
alpha_oa = 0.65*(exp((Y(i_V)+10)/-8.5)+exp((Y(i_V)-30)/-59.0))^(-1.0); %(Eq.37)
beta_oa = 0.65*(2.5+exp((Y(i_V)+82.0)/17.0))^(-1.0); %(Eq.37)
tau_oa = ((alpha_oa+beta_oa)^(-1.0))/K_Q10; %(Eq.38)
oa_infinity = (1.0+exp((Y(i_V)+20.47)/-17.54))^(-1.0); %(Eq.38)

% oi
alpha_oi = (18.53+1.0*exp((Y(i_V)+113.7)/10.95))^(-1.0); %(Eq.39)
beta_oi = (35.56+1.0*exp((Y(i_V)+1.26)/-7.44))^(-1.0); %(Eq.39)
tau_oi = ((alpha_oi+beta_oi)^(-1.0))/K_Q10; %(Eq.40)
oi_infinity = (1.0+exp((Y(i_V)+43.1)/5.3))^(-1.0); %(Eq.40)


%%% Ultrarapid Delayed Rectifier K+ Current
g_Kur = (0.005+0.05/(1.0+exp((Y(i_V)-15.0)/-13.0))); %(Eq.42)
IKur = Cm*fact_tot_IKur*g_Kur*Y(i_ua)^3.0*Y(i_ui)*(Y(i_V)-E_K); %(Eq.41)

% ua
alpha_ua = 0.65*(exp((Y(i_V)+10.0)/-8.5)+exp((Y(i_V)-30.0)/-59.0))^(-1.0); %(Eq.43)
beta_ua = 0.65*(2.5+exp((Y(i_V)+82.0)/17.0))^(-1.0); %(Eq.43)
tau_ua = ((alpha_ua+beta_ua)^(-1.0))/K_Q10; %(Eq.44)
ua_infinity = (1.0+exp((Y(i_V)+30.03)/-9.6))^(-1.0); %(Eq.44)

% ui
alpha_ui = (21.0+1.0*exp((Y(i_V)-185.0)/-28.0))^(-1.0); %(Eq.45)
beta_ui = exp((Y(i_V)-158.0)/(16.0));  %(Eq.45) same of 1.0/exp((Y(i_V)-158.0)/(-16.0));
tau_ui = (alpha_ui+beta_ui)^(-1.0)/K_Q10; %(Eq.46)
ui_infinity = (1.0+exp((Y(i_V)-99.45)/27.48))^(-1.0); %(Eq.46)

%%% Rapid Delayed Outward Rectifier K+ Current
IKr = Cm*fact_tot_IKr*g_Kr*Y(i_xr)*(Y(i_V)-E_K)/(1.0+exp((Y(i_V)+15.0)/22.4)); %(Eq.47)

% xr
if (abs(Y(i_V)+14.1) < 1.0e-10) %(Eq.48)
   alpha_xr = 0.0015;
else
   alpha_xr = 0.0003*(Y(i_V)+14.1)/(1.0-exp((Y(i_V)+14.1)/-5.0));
end
if (abs(Y(i_V)-3.3328) < 1.0e-10) %(Eq.48)
   beta_xr = 3.7836118e-4;
else
   beta_xr = 0.000073898*(Y(i_V)-3.3328)/(exp((Y(i_V)-3.3328)/5.1237)-1.0);
end
tau_xr = (alpha_xr+beta_xr)^(-1.0); %(Eq.49)
xr_infinity = (1.0+exp((Y(i_V)+14.1)/-6.5))^(-1.0); %(Eq.49)

%%% Slow Delayed Outward Rectifier K+ Current
IKs = Cm*fact_tot_IKs*g_Ks*Y(i_xs)^2.0*(Y(i_V)-E_K); %(Eq.50)

% xs
if (abs(Y(i_V)-19.9) < 1.0e-10) %(Eq.51)
   alpha_xs = 0.00068;
else
   alpha_xs = 0.00004*(Y(i_V)-19.9)/(1.0-exp((Y(i_V)-19.9)/-17.0));
end
if (abs(Y(i_V)-19.9) < 1.0e-10) %(Eq.51)
   beta_xs = 0.000315;
else
   beta_xs = 0.000035*(Y(i_V)-19.9)/(exp((Y(i_V)-19.9)/9.0)-1.0);
end
tau_xs = 0.5*(alpha_xs+beta_xs)^(-1.0); %(Eq.52)
xs_infinity = (1.0+exp((Y(i_V)-19.9)/-12.7))^(-0.5); %(Eq.52)

%%% I_KACh: Time-Independent, Acetylcholine-Activated K+ current (from Grandi 2011)
IKACh = Cm*fact_tot_IKACh*1/(1+(0.03/settings.ACh)^2.1) * (0.08 + 0.04/(1 + exp((Y(i_V)+91)/12))) * (Y(i_V) - E_K); % OK


%%% L-Type Ca2+ Current
ICaL = Cm*fact_tot_ICaL*g_Ca_L*Y(i_d)*Y(i_f)*Y(i_f_Ca)*(Y(i_V)-65.0); %(Eq.53)

% d
if (abs(Y(i_V)+10.0) < 1.0e-10) %(Eq.54)
   tau_d = 4.579/(1.0+exp((Y(i_V)+10.0)/-6.24)); 
else
   tau_d = (1.0-exp((Y(i_V)+10.0)/-6.24))/(0.035*(Y(i_V)+10.0)*(1.0+exp((Y(i_V)+10.0)/-6.24)));
end
d_infinity = (1.0+exp((Y(i_V)+10.0)/(-8.0)))^(-1.0); %(Eq.54)

% f
tau_f = 9.0*(0.0197*exp(-0.0337^2.0*(Y(i_V)+10.0)^2.0)+0.02)^(-1.0); %(Eq.55)
f_infinity = (1+exp((Y(i_V)+28)/6.9))^(-1); % (Eq.55) 

% fCa
tau_f_Ca = 2.0; %(Eq.56)
f_Ca_infinity = (1.0+Y(i_Ca_i)/0.00035)^(-1.0); %(Eq.56)

%%% Na+-K- Pump Current
sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0); %(Eq.59)
f_NaK = (1.0+0.1245*exp(-0.1*F*Y(i_V)/(R*T))+0.0365*sigma*exp(-F*Y(i_V)/(R*T)))^(-1.0); %(Eq.58)
INaK = Cm*fact_tot_INaK*i_NaK_max*f_NaK*(1.0/(1.0+(Km_Na_i/Y(i_Na_i))^1.5))*(K_o/(K_o+Km_K_o)); %(Eq.57)

%%% Na+/Ca2+ Exchanger Current
INaCa = Cm*fact_tot_INCX*I_NaCa_max*(exp(gamma*F*Y(i_V)/(R*T))*Y(i_Na_i)^3.0*Ca_o-exp((gamma-1.0)*F*Y(i_V)/(R*T))*Na_o^3.0*Y(i_Ca_i))/((K_mNa^3.0+Na_o^3.0)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*Y(i_V)*F/(R*T)))); %(Eq.60)

%%% Background Currents
Ib_Na = Cm*g_B_Na*(Y(i_V)-E_Na); %(Eq.62)
Ib_Ca = Cm*g_B_Ca*(Y(i_V)-E_Ca); %(Eq.61)
Ib_K  = Cm*g_B_K*(Y(i_V)-E_K); % null since g_B_K=0

%%% Ca2+ Pump Current
ICaP = Cm*fact_tot_ICaP*i_CaP_max*Y(i_Ca_i)/(0.0005+Y(i_Ca_i));

%%% Ca2+ Release Current from JSR
Irel = fact_tot_Irel*K_rel*Y(i_u)^2.0*Y(i_v)*Y(i_w)*(Y(i_Ca_rel)-Y(i_Ca_i)); %(Eq.64)
Fn = 1.0e3*(1.0e-15*V_rel*Irel-1.0e-15/(2.0*F)*(0.5*ICaL-0.2*INaCa));%(Eq.68)

% u
tau_u = 8.0; %(Eq.65)
u_infinity = 1/(1+exp(-(Fn-3.4175e-13)/13.67e-16));

% v
tau_v = 1.91+2.09*(1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^(-1.0); %(Eq.66)
v_infinity = 1.0-(1.0+exp(-(Fn-6.835e-14)/13.67e-16))^(-1.0); %(Eq.66)

% w
if (abs(Y(i_V)-7.9) < 1.0e-10) %(Eq.67)
   tau_w = 6.0*0.2/1.3; 
else
   tau_w = 6.0*(1.0-exp(-(Y(i_V)-7.9)/5.0))/((1.0+0.3*exp(-(Y(i_V)-7.9)/5.0))*1.0*(Y(i_V)-7.9));
end
w_infinity = 1.0-(1.0+exp(-(Y(i_V)-40.0)/17.0))^(-1.0); %(Eq.67)

%%% Transfer Current from NSR to JSR
Itr = fact_tot_Itr*(Y(i_Ca_up)-Y(i_Ca_rel))/tau_tr; %(Eq.69)

%%% Ca2+ Leak Current by the NSR
Iup_leak = fact_tot_Ileak*I_up_max*Y(i_Ca_up)/Ca_up_max; %(Eq.72)

%%% Ca2+ Uptake Current by the NSR
Iup = fact_tot_Iup*I_up_max/(1.0+K_up/Y(i_Ca_i)); %(Eq.71)

%%% Small Conductance Ca2+ - activated K+ channels: ISK
ISK = Cm*fact_tot_ISK*g_SK*Y(i_f_Ca_sk)^2.0*(Y(i_V)-E_K);

% fCa
tau_f_Ca_sk = 3;
q=2;
Ksk = 1.0000e-03;
f_Ca_infinity_sk = ((Y(i_Ca_i)/0.0025)^q)/(Ksk + ((Y(i_Ca_i)/0.0025)^q));

%%% K2P - VIOLETA 06/11/2024
IK2P = Cm*fact_tot_IK2P*g_K2P*Y(i_O_K2P)*(Y(i_V)-E_K);

% open gate for K2P channel
O_K2P_inf = 0.2 + (0.8/(1+exp(-(Y(i_V)-10)/14)));
O_K2P_tau = 2.0 + (40 / (1+exp((Y(i_V)+25)^2/50)));


%%% Stretch-activated current: ISAC
% I_SAC_Na
A_Na = P_Na*g_sac*Z_Na^2*F^2*Y(i_V)/(R*T);
B_Na = (Y(i_Na_i)-Na_o*exp(-Z_Na*F*Y(i_V)/(R*T)));
C_Na = (1-exp(-Z_Na*F*Y(i_V)/(R*T)));
ISAC_Na = fact_tot_ISAC * (A_Na*B_Na/C_Na);

% I_SAC_K
A_K = P_K*g_sac*Z_K^2*F^2*Y(i_V)/(R*T);
B_K = (Y(i_K_i)-K_o*exp(-Z_K*F*Y(i_V)/(R*T)));
C_K = (1-exp(-Z_K*F*Y(i_V)/(R*T)));
ISAC_K = fact_tot_ISAC * (A_K*B_K/C_K);

% I_SAC_Ca
A_Ca = P_Ca*g_sac*Z_Ca^2*F^2*Y(i_V)/(R*T);
B_Ca = (Y(i_Ca_i)-Ca_o*exp(-Z_Ca*F*Y(i_V)/(R*T)));
C_Ca = (1-exp(-Z_Ca*F*Y(i_V)/(R*T)));
ISAC_Ca = fact_tot_ISAC * (A_Ca*B_Ca/C_Ca);

% I_SAC
ISAC = ISAC_Na + ISAC_K + ISAC_Ca;

%%% Total ion current
Iion = INa + IK1 + Ito + IKur + IKr + IKs + Ib_Na + Ib_Ca + INaK ...
            + ICaP + INaCa + ICaL + IKACh + ISK + IK2P + ISAC;

ionCurrents = [INa, IK1, Ito, IKur, IKr, IKs, ICaL, INaK, INaCa, ICaP,...
                Ib_Na, Ib_Ca, Ib_K, Irel, Itr, Iup_leak, Iup, i_st, IKACh, ...
                ISK, IK2P, ISAC, ISAC_Na, ISAC_Ca, ISAC_K, Iion];
 

%-----------------------------------------------------------------------------
% Ca2+ Buffers
%-----------------------------------------------------------------------------

Ca_CMDN = CMDN_max*Y(i_Ca_i)/(Y(i_Ca_i)+Km_CMDN); %(Eq.73)
Ca_TRPN = TRPN_max*Y(i_Ca_i)/(Y(i_Ca_i)+Km_TRPN); %(Eq.74)
Ca_CSQN = CSQN_max*Y(i_Ca_rel)/(Y(i_Ca_rel)+Km_CSQN); %(Eq.75)


%-----------------------------------------------------------------------------
% Differential equations (update state variables)
%-----------------------------------------------------------------------------

dY(i_V) = -(Iion+i_st)/Cm;
dY(i_u)  = (u_infinity-Y(i_u))/tau_u;
dY(i_v)  = (v_infinity-Y(i_v))/tau_v;
dY(i_w)  = (w_infinity-Y(i_w))/tau_w;
dY(i_d)  = (d_infinity-Y(i_d))/tau_d;
dY(i_f_Ca)  = (f_Ca_infinity-Y(i_f_Ca))/tau_f_Ca;

dY(i_h)  = (h_inf-Y(i_h))/tau_h; 
dY(i_j)  = (j_inf-Y(i_j))/tau_j; 
dY(i_m)  = (m_inf-Y(i_m))/tau_m; 

dY(i_f_Ca_sk) = (f_Ca_infinity_sk -Y(i_f_Ca_sk))/tau_f_Ca_sk;
dY(i_O_K2P)   = (O_K2P_inf-Y(i_O_K2P))/O_K2P_tau; 

dY(i_f)  = (f_infinity-Y(i_f))/tau_f;
dY(i_xr) = (xr_infinity-Y(i_xr))/tau_xr;
dY(i_xs) = (xs_infinity-Y(i_xs))/tau_xs;
dY(i_oa) = (oa_infinity-Y(i_oa))/tau_oa;
dY(i_oi) = (oi_infinity-Y(i_oi))/tau_oi;
dY(i_ua) = (ua_infinity-Y(i_ua))/tau_ua;
dY(i_ui) = (ui_infinity-Y(i_ui))/tau_ui;

dY(i_Ca_up) = Iup-(Iup_leak+Itr*V_rel/V_up);
dY(i_Ca_rel) = (Itr-Irel)*(1.0+CSQN_max*Km_CSQN/(Y(i_Ca_rel)+Km_CSQN)^2.0)^(-1.0);

dY(i_Na_i) = (-3.0*INaK-(3.0*INaCa+Ib_Na+INa+ISAC_Na))/(V_i*F);
dY(i_K_i) = (2.0*INaK-(i_st+IK1+Ito+IKur+IKr+IKs+Ib_K+IKACh+ISK+IK2P+ISAC_K))/(V_i*F);

B1 = (2.0*INaCa-(ICaP+ICaL+Ib_Ca+ISAC_Ca))/(2.0*V_i*F)+(V_up*(Iup_leak-Iup)+Irel*V_rel)/V_i;
B2 = 1.0+TRPN_max*Km_TRPN/(Y(i_Ca_i)+Km_TRPN)^2.0+CMDN_max*Km_CMDN/(Y(i_Ca_i)+Km_CMDN)^2.0;
dY(i_Ca_i) = B1/B2; 


end
