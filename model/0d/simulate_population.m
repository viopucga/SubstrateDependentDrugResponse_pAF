function biomarkers = simulate_population(population_name, profile_ids, n_beats, bcl_ms)
%SIMULATE_POPULATION Simulate selected profiles and return biomarker tables.
if nargin < 3, n_beats = 50; end
if nargin < 4, bcl_ms = 1000; end
profile_ids = profile_ids(:);
rows = cell(numel(profile_ids),1);
for k = 1:numel(profile_ids)
    result = simulate_profile(population_name,profile_ids(k),n_beats,bcl_ms);
    ap = result.ap_biomarkers;
    ca = result.calcium_biomarkers;
    previous_t0 = max(0,(n_beats-2)*bcl_ms);
    [previous_t,previous_vm] = get_last(result.time_ms,result.vm_mV,previous_t0,(n_beats-1)*bcl_ms);
    previous_ap = get_AP_biomarkers(previous_t,previous_vm);
    alternans = 100*abs(previous_ap.APD90-ap.APD90)/mean([previous_ap.APD90,ap.APD90]) > 3;
    rows{k} = table(profile_ids(k),ap.RMP,ap.APA,ap.PlatA,ap.mx_dvdt,ap.APD95,ap.APD90,ap.APD50,ap.APD20,ap.rep_failure,ap.rep_EAD,ap.rep_DAD,alternans,ca.CaT_min,ca.CaT_max,ca.CaT_amp,ca.CaD90,ca.CaD50,ca.tau_decay,ca.SCaEs, ...
        'VariableNames',{'profile_id','RMP','APA','PlatA','mx_dvdt','APD95','APD90','APD50','APD20','rep_failure','rep_EAD','rep_DAD','alternans','CaT_min','CaT_max','CaT_amp','CaD90','CaD50','tau_decay','SCaEs'});
end
biomarkers = vertcat(rows{:});
end
