function population = GeneratePopulation_1(sim_config)
%GENERATEPOPULATION_1 Optionally generate a new Latin-hypercube population.
%
% Exact reproduction should use the supplied PoM_candidate_parameters.mat.
% Generation is disabled unless sim_config.allow_generation is true.

if ~isfield(sim_config,'allow_generation') || ~sim_config.allow_generation
    error('Candidate generation is disabled. Use the supplied population or explicitly set allow_generation=true.');
end
if ~isfield(sim_config,'size_pom'), sim_config.size_pom = 10000; end
if ~isfield(sim_config,'rng_seed')
    error('An explicit rng_seed is required when generating a new population.');
end

channel_names = {'GCaL','GK1','GKr','GKs','GKur','Gto','GKACh','GNa','GNaK','GNCX','GSK','GK2P','GSAC'};
rng(sim_config.rng_seed,'twister');
sample = lhsdesign(sim_config.size_pom,numel(channel_names));
low_limit = 0.5;
high_limit = 2.0;
parameter_factors = low_limit + (high_limit-low_limit).*sample;
population = array2table(parameter_factors,'VariableNames',channel_names);
end