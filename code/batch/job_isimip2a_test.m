% batch job for cluster: bsub -W 48:00 -R "rusage[mem=5000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
% MODULE:
%   isimip
% NAME:
%   job_isimip2a_test
% PURPOSE:
%   generate isimip flood damage with PAA=1 for all countries and save as a
%   csv file.
%   The job can be tested on a desktop, see run_on_desktop below and the
%   example.
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r ${climada_modules_folder}/isimip/code/batch/job_isimip2a_test.m bguillod@euler.ethz.ch:/cluster/home/bguillod/euler_jobs/
%   run on cluster:            bsub -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
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
% Benoit P. Guillod, bguillod@env.ethz.ch, 20180316, copy from
%   job_isimip_entities.m
%-


% PARAMETERS
% one should only have to edit this section
%
% if just a test with 2 small countries
TEST_ONLY=1;



run_on_desktop=strcmp(climada_machine_run,'local'); % to test the job on a desktop

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

% do one file for each model combination
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};

for iGHM=1:length(ghms)
    for iForcing=1:length(forcings)
        ghm=ghms{iGHM};
        forcing = forcings{iForcing};
    
        % check that this forcing-ghm combination exists
        [flddph_filename,~,fld_path] = isimip_get_flood_filename('2a', ghm, forcing, '0', 'historical');
        flood_filename=[fld_path filesep flddph_filename];
        if ~exist(flood_filename)
            continue
        end

        for i=1:length(all_countries)
            country=all_countries{i};
            output_i = isimip2a_FL_countrydata_for_PIK(country, ghm, forcing, params, params_damfun);
%             output = cat(2, all_years, repmat(string(country_iso3), [length(all_years) 1]),...
%                 repmat(string(continent), [length(all_years) 1]),...
%                 affected_area,  mean_flddph,...
%                 total_asset_value, total_asset_value_2005, ...
%                 exposed_asset_value, exposed_asset_value_2005, ...
%                 damage, damage_2005);

            if i==1
                output_all = output_i;
            else
                output_all = cat(1, output_all, output_i(2:end,:));
            end
            
        end
        output_file = [output_folder filesep 'output_' ghm '_' forcing '_FULL.csv'];
        output_all(ismissing(output_all))='NA';
        output_all2=cellstr(output_all);
        writetable(cell2table(output_all2),output_file,'writevariablenames',0);

        
    end
end

exit % the cluster appreciates this, gives back memory etc.
