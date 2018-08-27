% batch job for cluster: bsub -W 48:00 -R "rusage[mem=5000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_countryarea
% MODULE:
%   isimip
% NAME:
%   job_isimip2a_countryarea
% PURPOSE:
%   Determine country area of isimip assets and save as one file
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
%   bsub -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_countryarea
% EXAMPLE:
%   bsub -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_countryarea
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to work disk, see PARAMETERS
% MODIFICATION HISTORY:
% Benoit P. Guillod, bguillod@env.ethz.ch, 20180316, copy from
%   job_isimip_entities.m
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
else
    params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
    params.entity_prefix='FL1950';
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

for i=1:length(all_countries)
    country=all_countries{i};
    fprintf('\n   -> country: %s:\n', country);
    output_i = isimip2a_entities_area(country, params);
    if i==1
        output_all = [all_countries{i} output_i];
    else
        output_all = cat(1, output_all, [all_countries{i} output_i]);
    end
end
output_file = [output_folder filesep 'country_area_isimip.csv'];
output_all(ismissing(output_all))='NA';
output_all2=cellstr(output_all);
writetable(cell2table(output_all2),output_file,'writevariablenames',0);
fprintf('** Output file written: %s **\n', output_file);


exit % the cluster appreciates this, gives back memory etc.
