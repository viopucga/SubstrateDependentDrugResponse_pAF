%===============================================================================
%    Human atrial action potential model from Courtemanche et al. 1998 
%===============================================================================
%
% Courtemanche M, Ramirez RJ, Nattel S
% 1998 American Journal of Physiology-Heart and Circulatory Physiology, 275(1), H301-H321.
% PMID:9688927
% DOI: 10.1152/ajpheart.1998.275.1.H301
%
% Conversion from CellML 1.0 (www.cellml.org/) to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%
%-------------------------------------------------------------------------------
% Modifications to 1998 version (implemented and checked by Violeta Puche-García 2024):
%     - Factors to reproduce atrial electrical heterogeneity from Martinez-Mateu et al. 2018
%     - Factors used for the pAF and IREpAF populations in this study
%     - Fn variable solved as in cellml code and Sutanto et al. 2019 (MANTA), 
%       although it is not the equation presented in Courtemanche et al. 1998.
%     - ICaL variables (tau_d and f_infinity) solved as in cellml code and 
%       Courtemanche et al. 1998, although it is not the equation implemented in 
%       Sutanto et al. 2019 (MANTA)
%
% New currents added to the model
%     - Acetylcholine-activated K+ current (IKACh) adopted from Grandi et al. 2011
%       based on experimental data from Koumi et al. 1994 (added by Laura Martinez-Mateu 2018)
%            · Acetylcholine concentration is set through settings.ACh
%            · Acetylcholine concentration must be in [μM]
%     - SK current adopted from Celotto et al. 2023 (added by Violeta Puche-García 2024)
%            · SK channel conductance is set through settings.gSK
%            · For lone SK activity, gSK must be set to 0.0019 so APD90
%            variation is 12% from APD90 = 296ms for gSK=0 to APD90 = 251
%            for gK2P = 0.0091 (exp APD variation from Heijman net al. 2023)
%     - K2P current adopted from Schmidt et al. 2015 (added by Violeta Puche-García 2024)
%            · K2P channel conductance is set through settings.gK2P
%            · For lone K2P activity, gK2P must be set to 0.0091 so APD90
%            variation is 17% from APD90 = 296ms for gSK=0 to APD90 = 264
%            for gSK = 0.0019 (exp APD variation from Schmidt et al. 2015)
%     - Strech Activated Channel current (ISAC) adopted from Kuijpers et al. 2007 (added by Violeta Puche-García 2024)
%            · Channel conductance is set through settings.gSAC
%            · channel stretch is set through settings.SAC_lambda (for no strech SAC_lambda = 1)
%            · For lone SAC activity, gSAC must be set to 0.05 to reproduce
%            the results for increasing lamda in Kuijpers et al. 2007
%
% Simulation conditions (set through a struct named settings)
%     - Experiment configuration and stimulus parameters
%     - Drug study: concentration and IC50s-nH for ICaL, IK1, IKACh, IKr, IKs, IKur, Ito, ISK, IK2P, INa, INaK, INCX
%            · Concentration and IC50s must be in the same units
%     - Current factors (sensitivity, PoM): ICaL, IK1, IKACh, IKr, IKs, IKur, Ito, ISK, IK2P, INa, INaK, INCX, Irel, Ileak, Iup, Itr, ICaP
%
%       The default public configuration simulates the pAF left-atrial model.
%            · settings.AF_condition  = 'pAF';
%            · settings.atrial_model  = 'LA';
%            · settings.ACh           = 0.005;
%
%       Ioinc profile to simulate electrophysiological AP when combined activity
%       of IK,ACh (ACh = 0.005 microM), SK (gSK = 0.0019), K2P (gK2P = 0.0091),
%       SAC (gSAC = 0.05, SAC_lambda = 1.00)
%            · settings.fact_curr_ICaL  = 1.2;    · settings.fact_curr_INaK  = 1.0;
%            · settings.fact_curr_IKACh = 0.8;    · settings.fact_curr_INCX  = 1.2;
%            · settings.fact_curr_IK1   = 0.8;    · settings.fact_curr_Irel  = 0.7;
%            · settings.fact_curr_IKr   = 0.8;    · settings.fact_curr_Ileak = 0.7;
%            · settings.fact_curr_IKs   = 0.9;    · settings.fact_curr_Iup   = 1.2;
%            · settings.fact_curr_IKur  = 1.0;    · settings.fact_curr_Itr   = 0.7;
%            · settings.fact_curr_Ito   = 1.0;    · settings.fact_curr_ICaP  = 0.8;
%       Under this modifications, the following APD90 are obtained after 50p at BCL = 1000 ms:

%            · pAF-LA APD90 = 286 ms
%            · IREpAF-LA APD90 = 157 ms
%===============================================================================

function [t, vm, cai, StateVars, currents] = Courtemanche_main_all_PoM(settings,param_factors)

if nargin == 0                  % no initialized struct settings
    settings.marijose = 1;      % create a value to start the struct
end

settings = setDefaultSettings(settings);

y0       = getInitialVector(settings);

tol = 1e-12;
options  = odeset('AbsTol', tol,'RelTol', tol);

StateVars   = [];
Ti          = [];

t_start = tic;

if settings.comments == 1
    % verbose
    fprintf('Running simulation...  0%% completed');
end

% Get time (Ti) and state variables (StateVars)
for n = 1:settings.stim_num 

    if settings.comments == 1
        % verbose
        fprintf(repmat('\b', 1, length('100% completed')));
        fprintf('%3.0f%% completed', 100*(n/settings.stim_num))
    end

    % dY Integration
    [t Y] = ode15s('Courtemanche_model_all_PoM', [settings.stim_offset (settings.stim_BCL + settings.stim_offset)], y0, options, settings, param_factors);

    y0    = Y(end,:); % last row becomes initial vector for the following iteration

    Ti          = [Ti; t + (settings.stim_BCL*(n-1))];
    StateVars   = [StateVars; Y];

end

% Get currents only when requested by the caller
 currents = struct();
if nargout >= 5
Cm = 100; % pF
for i = 1:length(Ti)

        [dY, IsJs] = Courtemanche_model_all_PoM(Ti(i), StateVars(i,:), [], settings, param_factors); 

        currents.INa(i)         = IsJs(1)./Cm;     % [pA/pF]
        currents.IK1(i)         = IsJs(2)./Cm;     % [pA/pF]
        currents.Ito(i)         = IsJs(3)./Cm;     % [pA/pF]
        currents.IKur(i)        = IsJs(4)./Cm;     % [pA/pF]
        currents.IKr(i)         = IsJs(5)./Cm;     % [pA/pF]
        currents.IKs(i)         = IsJs(6)./Cm;     % [pA/pF]
        currents.ICa_L(i)       = IsJs(7)./Cm;     % [pA/pF]
        currents.INaK(i)        = IsJs(8)./Cm;     % [pA/pF]
        currents.INaCa(i)       = IsJs(9)./Cm;     % [pA/pF]
        currents.ICaP(i)        = IsJs(10)./Cm;    % [pA/pF]
        currents.IB_Na(i)       = IsJs(11)./Cm;    % [pA/pF]
        currents.IB_Ca(i)       = IsJs(12)./Cm;    % [pA/pF]
        currents.IB_K(i)        = IsJs(13)./Cm;    % [pA/pF]
        currents.Irel(i)        = IsJs(14)./Cm;    % [pA/pF]
        currents.Itr(i)         = IsJs(15)./Cm;    % [pA/pF]
        currents.Iup_leak(i)    = IsJs(16)./Cm;    % [pA/pF]
        currents.Iup(i)         = IsJs(17)./Cm;    % [pA/pF]
        currents.i_st(i)        = IsJs(18)./Cm;    % [pA/pF]
        currents.IKACh(i)       = IsJs(19)./Cm;    % [pA/pF]
        currents.ISK(i)         = IsJs(20)./Cm;    % [pA/pF]
        currents.IK2P(i)        = IsJs(21)./Cm;    % [pA/pF]
        currents.ISAC(i)        = IsJs(22)./Cm;    % [pA/pF]
        currents.ISAC_Na(i)     = IsJs(23)./Cm;    % [pA/pF]
        currents.ISAC_Ca(i)     = IsJs(24)./Cm;    % [pA/pF]
        currents.ISAC_K(i)      = IsJs(25)./Cm;    % [pA/pF]
        currents.Iion(i)        = IsJs(26)./Cm;    % [pA/pF]
        
end

end

% Re-organize output vectors
t = Ti';
vm = StateVars(:,15)';       % (mV) Membrane potential
cai = StateVars(:,10)'*1e+3; % (µM) Ca2+ intracellular concentration

t_elapsed = toc(t_start);
if settings.comments == 1
    fprintf('\n\t Simulation finished, took %.2f seconds \n', t_elapsed);
end



end



%%
function y0 = getInitialVector(settings)    
% Set the initial conditions for state variables

    y0 = [0	1	0.999000000000000	0.000137000000000000	0.775000000000000	0.999000000000000	0.965000000000000	0.978000000000000	0.00291000000000000	0.000102000000000000	1.49000000000000	1.49000000000000	139	11.2000000000000	-81.1800000000000	3.29000000000000e-05	0.0187000000000000	0.0304000000000000	0.999000000000000	0.00496000000000000	0.999000000000000 0.6403000000000000 0.2];

    if strcmp(settings.experiment,'patch_clamp') == 1
        % Vm(0) = -100
        y0 = [0	1	0.999000000000000	0.000137000000000000	0.775000000000000	0.999000000000000	0.965000000000000	0.978000000000000	0.00291000000000000	0.000102000000000000	1.49000000000000	1.49000000000000	139	11.2000000000000	-100.00000000000	3.29000000000000e-05	0.0187000000000000	0.0304000000000000	0.999000000000000	0.00496000000000000	0.999000000000000 0.6403000000000000 0.2];
    end

end



%%
function settings = setDefaultSettings(settings)
% Set default settings for the simulation
   
    if ~isfield(settings, 'comments'),              settings.comments               = 1;                   end      % 0 if no comments, 1 if so

    % Stimulus settings from the paper [ Ref. CTR (1998) ]
	if ~isfield(settings, 'stim_num'),              settings.stim_num               = 12;                 end      % # of stimuli
	if ~isfield(settings, 'stim_offset'),           settings.stim_offset            = 0.0;                end	   % offset for the stimuls (ms)
    if ~isfield(settings, 'stim_amp'),              settings.stim_amp               = -2000;              end      % Stimulation amplitude (pA/pF)
    if ~isfield(settings, 'stim_dur'),              settings.stim_dur               = 2;                  end      % Stimulation duration (ms)
    if ~isfield(settings, 'stim_BCL'),              settings.stim_BCL               = 1000;               end      % Basic Cycle Length (ms)
    if ~isfield(settings, 'Tsim'),                  settings.Tsim                   = 12000;              end      % Duration of simulation (ms)
           
    % Experiment configuration
    if ~isfield(settings, 'experiment'),            settings.experiment             = 'AP';               end      % experiment [AP, patch_clamp]
    if ~isfield(settings, 'AF_condition'),          settings.AF_condition           = 'pAF';              end      % Study population [pAF, IREpAF]
    if ~isfield(settings, 'atrial_model'),          settings.atrial_model           = 'LA';            end      % atrial model ['RA_PM', 'CT_BBra', 'TVR', 'RAA', 'LA', 'LAA', 'PV', 'MVR', 'BBla']
    if ~isfield(settings, 'ACh'),                   settings.ACh                    = 0.0;                end      % ACh concentration [µM]
    if ~isfield(settings, 'gSK'),                   settings.gSK                    = 0.0;                end      % conductivity of SK channel
    if ~isfield(settings, 'gK2P'),                  settings.gK2P                   = 0.0;                end      % conductivity of K2P channel
    if ~isfield(settings, 'gSAC'),                  settings.gSAC                   = 0.0;                end      % maximun membrane conductance
    if ~isfield(settings, 'SAC_lambda'),            settings.SAC_lambda             = 0.0;                end      % strech ratio (1 if cell is not stretched)

    % Drug parameters (µM)
    if ~isfield(settings, 'drug_conc'),             settings.drug_conc               = 0;                 end      % Drug concentration (µM)
    if ~isfield(settings, 'drug_ICaL_IC50'),        settings.drug_ICaL_IC50          = 1e+6;              end      % ICaL_IC50 (µM)
    if ~isfield(settings, 'drug_ICaL_nH'),          settings.drug_ICaL_nH            = 1;                 end           
    if ~isfield(settings, 'drug_IK1_IC50'),         settings.drug_IK1_IC50           = 1e+6;              end      % IK1_IC50 (µM)
    if ~isfield(settings, 'drug_IK1_nH'),           settings.drug_IK1_nH             = 1;                 end           
    if ~isfield(settings, 'drug_IKACh_IC50'),       settings.drug_IKACh_IC50         = 1e+6;              end      % IKACh_IC50 (µM)
    if ~isfield(settings, 'drug_IKACh_nH'),         settings.drug_IKACh_nH           = 1;                 end 
    if ~isfield(settings, 'drug_IKr_IC50'),         settings.drug_IKr_IC50           = 1e+6;              end      % IKr_IC50 (µM)
    if ~isfield(settings, 'drug_IKr_nH'),           settings.drug_IKr_nH             = 1;                 end        
    if ~isfield(settings, 'drug_IKs_IC50'),         settings.drug_IKs_IC50           = 1e+6;              end      % IKs_IC50 (µM)
    if ~isfield(settings, 'drug_IKs_nH'),           settings.drug_IKs_nH             = 1;                 end          
    if ~isfield(settings, 'drug_IKur_IC50'),        settings.drug_IKur_IC50          = 1e+6;              end      % IKur_IC50 (µM)
    if ~isfield(settings, 'drug_IKur_nH'),          settings.drug_IKur_nH            = 1;                 end      
    if ~isfield(settings, 'drug_Ito_IC50'),         settings.drug_Ito_IC50           = 1e+6;              end      % Ito_IC50 (µM)
    if ~isfield(settings, 'drug_Ito_nH'),           settings.drug_Ito_nH             = 1;                 end      
    if ~isfield(settings, 'drug_ISK_IC50'),         settings.drug_ISK_IC50           = 1e+6;              end      % ISK_IC50 (µM)
    if ~isfield(settings, 'drug_ISK_nH'),           settings.drug_ISK_nH             = 1;                 end 
    if ~isfield(settings, 'drug_IK2P_IC50'),        settings.drug_IK2P_IC50          = 1e+6;              end      % IK2P_IC50 (µM)
    if ~isfield(settings, 'drug_IK2P_nH'),          settings.drug_IK2P_nH            = 1;                 end 
    if ~isfield(settings, 'drug_INa_IC50'),         settings.drug_INa_IC50           = 1e+6;              end      % INa_IC50 (µM)
    if ~isfield(settings, 'drug_INa_nH'),           settings.drug_INa_nH             = 1;                 end      
    if ~isfield(settings, 'drug_INaK_IC50'),        settings.drug_INaK_IC50          = 1e+6;              end      % INaK_IC50 (µM)
    if ~isfield(settings, 'drug_INaK_nH'),          settings.drug_INaK_nH            = 1;                 end      
    if ~isfield(settings, 'drug_INCX_IC50'),        settings.drug_INCX_IC50          = 1e+6;              end      % INCX_IC50 (µM)
    if ~isfield(settings, 'drug_INCX_nH'),          settings.drug_INCX_nH            = 1;                 end      
         
    % Current factors
    if ~isfield(settings, 'fact_curr_ICaL'),        settings.fact_curr_ICaL          = 1;                 end      % factor to ICaL conductance  
    if ~isfield(settings, 'fact_curr_IK1'),         settings.fact_curr_IK1           = 1;                 end      % factor to IK1 conductance   
    if ~isfield(settings, 'fact_curr_IKACh'),       settings.fact_curr_IKACh         = 1;                 end      % factor to IKACh conductance 
    if ~isfield(settings, 'fact_curr_IKr'),         settings.fact_curr_IKr           = 1;                 end      % factor to IKr conductance   
    if ~isfield(settings, 'fact_curr_IKs'),         settings.fact_curr_IKs           = 1;                 end      % factor to IKs conductance   
    if ~isfield(settings, 'fact_curr_IKur'),        settings.fact_curr_IKur          = 1;                 end      % factor to IKur conductance  
    if ~isfield(settings, 'fact_curr_Ito'),         settings.fact_curr_Ito           = 1;                 end      % factor to Ito conductance   
    if ~isfield(settings, 'fact_curr_ISK'),         settings.fact_curr_ISK           = 1;                 end      % factor to ISK conductance   
    if ~isfield(settings, 'fact_curr_IK2P'),        settings.fact_curr_IK2P          = 1;                 end      % factor to IK2P conductance
    if ~isfield(settings, 'fact_curr_ISAC'),        settings.fact_curr_ISAC          = 1;                 end      % factor to ISAC conductance
    if ~isfield(settings, 'fact_curr_INa'),         settings.fact_curr_INa           = 1;                 end      % factor to INa conductance   
    if ~isfield(settings, 'fact_curr_INaK'),        settings.fact_curr_INaK          = 1;                 end      % factor to INaK conductance  
    if ~isfield(settings, 'fact_curr_INCX'),        settings.fact_curr_INCX          = 1;                 end      % factor to INCX conductance  
    if ~isfield(settings, 'fact_curr_Irel'),        settings.fact_curr_Irel          = 1;                 end      % factor to Irel conductance  
    if ~isfield(settings, 'fact_curr_Ileak'),       settings.fact_curr_Ileak         = 1;                 end      % factor to Ileak conductance 
    if ~isfield(settings, 'fact_curr_Iup'),         settings.fact_curr_Iup           = 1;                 end      % factor to Iup conductance   
    if ~isfield(settings, 'fact_curr_Itr'),         settings.fact_curr_Itr           = 1;                 end      % factor to Itr conductance   
    if ~isfield(settings, 'fact_curr_ICaP'),        settings.fact_curr_ICaP          = 1;                 end      % factor to ICaP conductance  

end


