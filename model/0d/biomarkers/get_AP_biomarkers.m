function biomarkers = get_AP_biomarkers(last_AP_t,last_AP_Vm)

%====================================================================================
% get_AP_biomarkers
%       Outputs:
%           RMP (mV)        Resting Membrane Potential
%           APO (mV)    	from 0 to Vmax
%           APA (mV)	    amplitude from RMP to Vmax
%           mx_dvdt (V/s)   maximun depolarization velocity
%           APD90 (ms)      AP duration at 90% repolarization
%           APD50 (ms)      AP duration at 50% repolarization
%
%       Inputs:
%           t (ms)          time instants for the last AP
%           Vm (mV)         membrane potential for the last AP
%
%====================================================================================

% ------ WARNINGS ---------

if isempty(last_AP_t) || isempty(last_AP_Vm)
    warning('t and Vm should not be empty');
    return;
end

if length(last_AP_t) ~= length(last_AP_Vm)
    warning('t and Vm must be the same length');
    return;
end


% ------ BIOMARKERS ---------

% RMP
RMP = last_AP_Vm(1);
if last_AP_Vm(end) < last_AP_Vm(1)
    RMP         = last_AP_Vm(end);
end

% APO
[max_Vm,idx_max_Vm] = findpeaks(last_AP_Vm,'MinPeakProminence',10);
max_Vm=max_Vm(1);
idx_max_Vm=idx_max_Vm(1);
APO = max_Vm(1);

% APA
APA = APO - RMP;

% mx_dvdt
dv              = diff(last_AP_Vm);
dt              = diff(last_AP_t);
dvdt            = dv./dt;
dvdt(dvdt==Inf)=0;
[mx_dvdt, idx_mx_dvdt] = max(dvdt);

% APD90
V95 = APO - (APA*0.95);
APD95_f = idx_max_Vm;
while last_AP_Vm(APD95_f) > V95
    
    if APD95_f == length(last_AP_Vm), break % when APD90 matched BCL
    else
        APD95_f = APD95_f + 1;
    end
    
end
APD95 = last_AP_t(APD95_f) - last_AP_t(idx_mx_dvdt);

% APD90
V90 = APO - (APA*0.9);
APD90_f = idx_max_Vm;
while last_AP_Vm(APD90_f) > V90
    
    if APD90_f == length(last_AP_Vm), break % when APD90 matched BCL
    else
        APD90_f = APD90_f + 1;
    end
    
end
APD90 = last_AP_t(APD90_f) - last_AP_t(idx_mx_dvdt);

% APD50
V50 = APO - (APA*0.5);
APD50_f = idx_max_Vm;
while last_AP_Vm(APD50_f) > V50
    if APD50_f == length(last_AP_Vm), break % when APD50 matched BCL
    else
        APD50_f = APD50_f + 1;
    end
end
APD50 = last_AP_t(APD50_f) - last_AP_t(idx_mx_dvdt);

% APD_30
V30 = APO - (APA*0.3);
APD30_f = idx_max_Vm;
while last_AP_Vm(APD30_f) > V30
    if APD30_f == length(last_AP_Vm), break % when APD30 matched BCL
    else
        APD30_f = APD30_f + 1;
    end
end
APD30 = last_AP_t(APD30_f) - last_AP_t(idx_mx_dvdt);

% APD_20
V20 = APO - (APA*0.2);
APD20_f = idx_max_Vm;
while last_AP_Vm(APD20_f) > V20
    if APD20_f == length(last_AP_Vm), break % when APD20 matched BCL
    else
        APD20_f = APD20_f + 1;
    end
end
APD20 = last_AP_t(APD20_f) - last_AP_t(idx_mx_dvdt);

% PLT_20 (Wetter, 2013 - nSR/cAF human RAA recordings) 
% plateau potential defined as the mean potential (mV) in the time window between 20 and 30% of repolarization
PLT_20 = mean(last_AP_Vm(APD20_f:APD30_f));

% PlatA (Verkerk 2021 - human LAA recordings)
% AP plateau amplitude measured 20 ms after initiation of the AP upstroke
idx_tPlatA = find(last_AP_t>last_AP_t(idx_max_Vm)+20); idx_tPlatA = idx_tPlatA(1);  % idx for t= Vmax + 20ms
PlatA = last_AP_Vm(idx_tPlatA) - RMP;

% repolarization abnormalities
% EAD - at least 2 local maxima between APO and APD90 
%       (jllopis: any event with positive voltage gradient after 100ms from the beggining of AP
%        =/ pueden darse AP con pico de joroba posterior a los 100ms)
% DAD - any event with positive voltage gradient after APD90 from the beggining of AP
% repolarization failure - Vm at end of simulation > -40mV
% rep_EAD_Vm = last_AP_Vm(idx_max_Vm:APD90_f);
TF = find_local_max(last_AP_Vm(idx_max_Vm:APD90_f)); % [TF,~] = islocalmax(last_AP_Vm(idx_max_Vm:APD90_f)); %-> no funciona en BSC por ser MATLAB 2024
rep_EAD = numel(find(TF==1)) > 1; % 0 if no EAD, 1 if there is EAD
% idx_t100 = find(last_AP_t>last_AP_t(1)+100); idx_t100 = idx_t100(1);  % idx for t=100ms
% rep_EAD = sum(unique(dvdt(idx_t100:end) > 0)); % 0 if no EAD, 1 if there is EAD
rep_DAD = numel(unique(dvdt(APD90_f:end) > 0)) > 1; % 0 if no DAD, 1 if there is DAD
rep_failure = (last_AP_Vm(end) > -40);


% 2struct
biomarkers.RMP            = RMP;
% biomarkers.APO            = APO;
biomarkers.APA            = APA;
biomarkers.PlatA          = PlatA;
biomarkers.mx_dvdt        = mx_dvdt;
biomarkers.APD95          = APD95;
biomarkers.APD90          = APD90;
biomarkers.APD50          = APD50;
% biomarkers.APD30          = APD30;
biomarkers.APD20          = APD20;
% biomarkers.triang         = APD90 - APD30;
% biomarkers.PLT_20         = PLT_20;
% biomarkers.APD90_f        = APD90_f;
% biomarkers.APD90_i        = idx_mx_dvdt;
% biomarkers.idx_max_Vm     = idx_max_Vm;
biomarkers.rep_failure    = rep_failure;
biomarkers.rep_EAD        = rep_EAD;
biomarkers.rep_DAD        = rep_DAD;

end

