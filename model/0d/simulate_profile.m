function result = simulate_profile(population_name, profile_id, n_beats, bcl_ms)
%SIMULATE_PROFILE Simulate one candidate profile with the study 0D model.
%
% result = simulate_profile(population_name, profile_id)
% result = simulate_profile(population_name, profile_id, n_beats, bcl_ms)
%
% population_name must be "pAF" or "IREpAF". profile_id is the stable ID
% used in data/populations/population_10000.csv. The default protocol uses
% 50 beats at a basic cycle length of 1000 ms.

if nargin < 3, n_beats = 50; end
if nargin < 4, bcl_ms = 1000; end
population_name = string(population_name);
if ~ismember(population_name, ["pAF", "IREpAF"])
    error('population_name must be pAF or IREpAF.');
end

model_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(model_dir));
addpath(model_dir, fullfile(model_dir, 'biomarkers'));

population_file = fullfile(repo_root, 'data', 'populations', 'PoM_candidate_parameters.mat');
loaded = load(population_file);
population = loaded.PoM_candidate_parameters;
if profile_id < 1 || profile_id > height(population) || profile_id ~= fix(profile_id)
    error('profile_id must be an integer between 1 and %d.', height(population));
end

conductance_names = {'GCaL','GK1','GKr','GKs','GKur','Gto','GKACh','GNa','GNaK','GNCX','GSK','GK2P','GSAC'};
parameter_factors = population{profile_id, conductance_names};
baseline_factors = [1.2, 0.8, 0.8, 0.9, 1.0, 1.0, 0.8, 1.0, 1.0, 1.2, 1.0, 1.0, 1.0];
parameter_factors = parameter_factors .* baseline_factors;

settings.comments = 0;
settings.stim_num = n_beats;
settings.stim_BCL = bcl_ms;
settings.AF_condition = char(population_name);
settings.atrial_model = 'LA';
settings.ACh = 0.005;
settings.gSK = 0.0019;
settings.gK2P = 0.0091;
settings.gSAC = 0.05;
settings.SAC_lambda = 1.0;
settings.fact_curr_Irel = 0.7;
settings.fact_curr_Ileak = 0.7;
settings.fact_curr_Iup = 1.2;
settings.fact_curr_Itr = 0.7;
settings.fact_curr_ICaP = 0.8;

[t, vm, cai] = Courtemanche_main_all_PoM(settings, parameter_factors);
[last_t, last_vm] = get_last(t, vm, (n_beats-1)*bcl_ms, n_beats*bcl_ms);
[~, last_cai] = get_last(t, cai, (n_beats-1)*bcl_ms, n_beats*bcl_ms);

result.population = population_name;
result.profile_id = profile_id;
result.n_beats = n_beats;
result.bcl_ms = bcl_ms;
result.parameter_factors = parameter_factors;
result.time_ms = t;
result.vm_mV = vm;
result.cai_uM = cai;
result.last_beat.time_ms = last_t;
result.last_beat.vm_mV = last_vm;
result.last_beat.cai_uM = last_cai;
result.ap_biomarkers = get_AP_biomarkers(last_t, last_vm);
result.calcium_biomarkers = get_CaT_biomarkers(last_t, last_cai);
end
