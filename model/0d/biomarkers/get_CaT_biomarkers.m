function biomarkers = get_CaT_biomarkers(last_AP_t,last_AP_cai)

%====================================================================================
% get_AP_biomarkers
%       Outputs:
%           max CaT (µM)        Maxium CaT peak
%           resting CaT (µM)    Resting CaT concentration
%           amplitude (µM)	    amplitude from resting to CaTmax values
%           CaD50 (ms)          Time to reach 50% recovery
%           CaD90 (ms)          Time to reach 90% recovery
%           Slope ?????
%
%       Inputs:
%           t (ms)              time instants for the last AP
%           Cai (µM)            intracellular calcium transient for the last AP
%
%====================================================================================

% ------ WARNINGS ---------

if isempty(last_AP_t) || isempty(last_AP_cai)
    warning('t and cai should not be empty');
    return;
end

if length(last_AP_t) ~= length(last_AP_cai)
    warning('t and cai must be the same length');
    return;
end


% ------ BIOMARKERS ---------

% resting value
CaT_min = last_AP_cai(1);
if last_AP_cai(end) < last_AP_cai(1)
    CaT_min         = last_AP_cai(end);
end

% maximum peak
[CaT_max,idx_max_CaT] = max(last_AP_cai);

% amplitude
CaT_amp = CaT_max - CaT_min;

% CaD90
cat90 = CaT_max - (CaT_amp*0.9);
CaD90_f = idx_max_CaT;
while last_AP_cai(CaD90_f) > cat90
    
    if CaD90_f == length(last_AP_cai), break % when APD90 matched BCL
    else
        CaD90_f = CaD90_f + 1;
    end
    
end
CaD90 = last_AP_t(CaD90_f) - last_AP_t(idx_max_CaT);

% CaD50
cat50 = CaT_max - (CaT_amp*0.5);
CaD50_f = idx_max_CaT;
while last_AP_cai(CaD50_f) > cat50
    if CaD50_f == length(last_AP_cai), break % when APD50 matched BCL
    else
        CaD50_f = CaD50_f + 1;
    end
end
CaD50 = last_AP_t(CaD50_f) - last_AP_t(idx_max_CaT);

% tau decay
time = last_AP_t-last_AP_t(1);
ca_i = last_AP_cai; 

dcdt = gradient(ca_i, time);        % first derivative
d2cdt2 = gradient(dcdt, time);      % second derivative


% idx_t50 = find(time>50); idx_t50 = idx_t50(1);  % idx for t=50
% [~,inflex_point] = max(d2cdt2(idx_t50:end));    % point where second derivative cuts zero
% inflex_point= idx_t50+inflex_point;
% t_decay = time(inflex_point:end);     % Select the decay phase data from the concave/convexe to the end
% ca_i_decay = ca_i(inflex_point:end);
% 
% [~, peak_index] = max(ca_i);
% t_decay = time(peak_index:end);
% ca_i_decay = ca_i(peak_index:end);
% 
% idx_t300 = find(time>100); idx_t300 = idx_t300(1);  % idx for t=50
% t_decay = time(idx_t300:end);     % Select the decay phase data from the concave/convexe to the end
% ca_i_decay = ca_i(idx_t300:end);
% 
% t_decay = time(CaD50_f:end);     % Select the decay phase data from the concave/convexe to the end
% ca_i_decay = ca_i(CaD50_f:end);

[~,inflex_point] = max(d2cdt2(idx_max_CaT:end));
inflex_point = inflex_point + idx_max_CaT;
t_decay = time(inflex_point:end);     % Select the decay phase data from the concave/convexe to the end
ca_i_decay = ca_i(inflex_point:end);

exp_decay = @(params, t) params(1) * exp(-t / params(2)) + params(3);       % Define the exponential decay function
initial_guess = [CaT_amp, 250, CaT_min];                                    % Initial guess for the parameters [CaT_amp, tau, CaT_diast]
obj_fun = @(params) sum((ca_i_decay - exp_decay(params, t_decay)).^2);      % Define the objective function for fitting

options = optimset('Display', 'off');                                       % Use fminsearch to fit the exponential decay model (suppress output display)
fitted_params = fminsearch(obj_fun, initial_guess, options);
tau_decay = fitted_params(2);

% % Plot the fitted results
% fitted_curve = exp_decay(fitted_params, t_decay);
% 
% figure; 
% subplot(2,1,1); hold on;
% plot(time, ca_i, 'k-', 'DisplayName', 'CaT evolution');
% plot(t_decay, ca_i_decay, 'b-', 'DisplayName', 'Data to fit');
% plot(t_decay, fitted_curve, 'r--', 'DisplayName', ['Fitted Curve (tau_{decay} = ' num2str(round(tau_decay,1)) ' ms)']);
% ylabel('[Ca^{2+}]_i');
% legend boxoff;
% 
% subplot(2,1,2); hold on;
% plot(time, ca_i, 'b', 'DisplayName', '[Ca^{2+}]_i');
% plot(time, d2cdt2, 'g', 'DisplayName', 'Segunda derivada d^2[Ca^{2+}]_i/dt^2');
% plot(time(inflex_point),ca_i(inflex_point),'*r','DisplayName','punto inflexion')
% ylim([-0.1 1]); xlabel('time (ms)')
% legend;
% set(gcf,'Position',[2059         240         560         643])

% repolarization abnormalities
% EAD - any event with positive voltage gradiente after 100ms from the beggining of AP
idx_t300 = find(last_AP_t>last_AP_t(1)+300); idx_t300 = idx_t300(1);  % idx for t=300ms
SCaEs = sum(unique(dcdt(idx_t300:end) > 0)); % 0 if no EAD, 1 if there is EAD

% 2struct
biomarkers.CaT_min        = CaT_min;
biomarkers.CaT_max        = CaT_max;
biomarkers.CaT_amp        = CaT_amp;
biomarkers.CaD90          = CaD90;
biomarkers.CaD50          = CaD50;
biomarkers.tau_decay      = tau_decay;
biomarkers.SCaEs          = SCaEs;


end

