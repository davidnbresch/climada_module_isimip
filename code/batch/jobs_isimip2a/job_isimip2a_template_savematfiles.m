% batch job for cluster: bsub -W 48:00 -R "rusage[mem=5000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_template
% MODULE:
%   isimip
% NAME:
%   job_isimip2a_test
% PURPOSE:
%   PROPER DOCUMENTATION NEEDS UPDATING HERE
%   generate isimip flood damage with PAA=1 for all countries and save as a
%   csv file.
%   The job can be tested on a desktop, see run_on_desktop below and the
%   example.
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r ${climada_modules_folder}/isimip/code/batch/jobs_isimip2a/job_isimip2a_template.m bguillod@euler.ethz.ch:/cluster/home/bguillod/euler_jobs/
%   copy script to launch all jobs:  scp -r ${climada_modules_folder}/isimip/code/batch/jobs_isimip2a/LAUNCH_ALL.sh bguillod@euler.ethz.ch:/cluster/home/bguillod/euler_jobs/

%   run on cluster (mem=20000, so it can work for Canada, Russia etc):
%                              bsub -R "rusage[mem=20000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
%   or run all of them using ./LAUNCH_ALL.sh
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/entities/*.mat Documents/_GIT/climada_data/entities/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
% EXAMPLE:
%   bsub -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to work disk, see PARAMETERS
% MODIFICATION HISTORY:
% Benoit P. Guillod, bguillod@env.ethz.ch, 20180718, initial based on
% job_isimip2a_template.m
%-


% PARAMETERS
% one should only have to edit this section
%
% if just a test with 2 small countries
TEST_ONLY=0;

% are we on euler?
if strcmpi('GLNXA64',computer) % on Euler, 20180101
    run_on_desktop=0;
else
    run_on_desktop=1;
end

if ~run_on_desktop
    cd /cluster/home/bguillod/climada % to make sure the cluster finds climada
    startup % climada_global exists afterwards
    pwd % just to check where the job is running from
    % N_pool_workers=24; % for parpool
    % climada_global.parfor=1; % for parpool
end


% pool=parpool(N_pool_workers);

% define where data will be saved
output_folder=[climada_global.data_dir filesep 'isimip/results'];

clear params;
% parameters on cluster
if run_on_desktop
    % local parameters (for testing only)
%     params.entity_folder=climada_global.entities_dir;
    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];
    params.entity_prefix='FL1950';
    params_damfun.filename_suffix='PAA1';
    params_damfun.filepath='/Users/bguillod/Documents/work/ETH/floods/damage_functions/files_mine';
else
    params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
    params.entity_prefix='FL1950';
    params_damfun.filename_suffix='PAA1';
    params_damfun.filepath=[climada_global.data_dir filesep 'isimip/entities/damfun'];
end

% generate a list of all countries to loop over
NatID_RegID=isimip_NatID_RegID; % get the mapping ISO3 - isimip country code
all_countries = NatID_RegID.ISO3;
% test: run only a subset of countries
if TEST_ONLY
    %     all_countries = all_countries([1 5]);
    all_countries = all_countries(1:10);
end

% keep only countries which 'exist' as ISO3
all_countries_keep=ones(1,length(all_countries));
for i=1:length(all_countries)
    country=all_countries{i};
    country_exists =  climada_country_name(country);
    if isempty(country_exists)
        all_countries_keep(i)=0;
    end
end
fprintf('WARNING: Skipping the following %s countries:\n', num2str(sum(~all_countries_keep)));
all_countries(~all_countries_keep)
all_countries=all_countries(find(all_countries_keep));
clear all_countries_keep country_exists;

ghm='GHMNAME';
forcing = 'FORCINGNAME';
subtract_matsiro = MATSIROSUBTRACTION;
if subtract_matsiro
    append_name_matsiro='_mFRCmastiro';
else
    append_name_matsiro='';
end

fprintf('\n***** Running model %s-%s *****\n', ghm, forcing);

% check that this forcing-ghm combination exists
[flddph_filename,~,fld_path] = isimip_get_flood_filename('2a', ghm, forcing, '0', 'historical');
flood_filename=[fld_path filesep flddph_filename];
if ~exist(flood_filename)
    %          continue
    fprintf('\n   ** THIS MODEL COMBINATION DOES NOT EXIST, CANCEL SCRIPT: %s:\n', flood_filename);
    return
end

for i=1:length(all_countries)
    country=all_countries{i};
    fprintf('\n   -> country: %s:\n', country);
    % temporary output folder for individual country data
    output_folder_temp=[output_folder filesep 'step_data' filesep ghm '_' forcing];
    if ~exist(output_folder_temp, 'dir')
        mkdir(output_folder_temp);
    end
    file_out = [output_folder_temp filesep ghm '_' forcing append_name_matsiro '_' country '.mat'];
    % if temporary file for this country exists, go to next
    if exist(file_out, 'file')
        continue
    end
    output_i = isimip2a_FL_countrydata_for_PIK(country, ghm, forcing, params, params_damfun, subtract_matsiro);
    % save temporary file with country data
    save(file_out,'output_i');
    fprintf('** Output file for one country written: %s **\n', file_out);
end


exit % the cluster appreciates this, gives back memory etc.
