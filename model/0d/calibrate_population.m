function [accepted, checks] = calibrate_population(biomarkers_bcl1000, biomarkers_bcl500, ranges_file)
%CALIBRATE_POPULATION Apply the study calibration criteria to 0D biomarkers.
%
% [accepted, checks] = calibrate_population(b1000, b500)
%
% b1000 and b500 may be tables or MAT filenames containing a table named
% `bmkrs`. The default ranges are read from data/populations/bmkrs_paper.mat.

model_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(model_dir));
addpath(fullfile(model_dir, 'biomarkers'));
if nargin < 3
    ranges_file = fullfile(repo_root, 'data', 'populations', 'bmkrs_paper.mat');
end

b1000 = load_biomarker_table(biomarkers_bcl1000);
b500 = load_biomarker_table(biomarkers_bcl500);
if height(b1000) ~= height(b500)
    error('BCL 1000 and BCL 500 tables must contain the same profiles.');
end

ranges_data = load(ranges_file);
ranges = ranges_data.bmkrs_ref;
range_of = @(name) [ranges{strcmp(ranges(:,1),name),2}, ranges{strcmp(ranges(:,1),name),3}];

checks = table();
checks.APD90 = in_range(b1000.APD90, range_of('APD90'));
checks.APD50 = in_range(b1000.APD50, range_of('APD50'));
checks.APD20 = in_range(b1000.APD20, range_of('APD20'));
checks.max_upstroke = in_range(b1000.mx_dvdt, range_of('mx_dvdt'));
checks.plateau_amplitude = in_range(b1000.PlatA, range_of('PlatA'));
checks.APA = in_range(b1000.APA, range_of('APA'));
checks.RMP = in_range(b1000.RMP, range_of('RMP'));
checks.CaT_min = in_range(b500.CaT_min, range_of('CaT_min'));
checks.CaT_max = in_range(b500.CaT_max, range_of('CaT_max'));
checks.CaT_amplitude = in_range(b500.CaT_amp, range_of('CaT_amp'));
checks.tau_decay = in_range(b500.tau_decay, range_of('tau_decay'));
checks.no_repolarization_failure = ~logical(b500.rep_failure);
checks.no_EAD = ~logical(b500.rep_EAD);
accepted = all(checks{:,:},2);
checks.accepted = accepted;
end

function biomarkers = load_biomarker_table(input_value)
if istable(input_value)
    biomarkers = input_value;
elseif ischar(input_value) || isstring(input_value)
    data = load(input_value);
    if ~isfield(data,'bmkrs') || ~istable(data.bmkrs)
        error('MAT file must contain a biomarker table named bmkrs.');
    end
    biomarkers = data.bmkrs;
else
    error('Biomarkers must be supplied as a table or MAT filename.');
end
end
