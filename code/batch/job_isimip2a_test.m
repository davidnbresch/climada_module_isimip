% batch job for cluster: bsub -W 48:00 -R "rusage[mem=5000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2a_test
% MODULE:
%   isimip
% NAME:
%   job_isimip2a_test
% PURPOSE:
%   generate isimip flood damage with PAA=1 for all countries and save as a
%   csv file.
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/batch/job_isimip_entities.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   run on cluster:            bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip_entities
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/entities/*.mat Documents/_GIT/climada_data/entities/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_ispwdimip_entities
% EXAMPLE:
%   bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip_entities
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
cd /cluster/home/bguillod/climada % to make sure the cluster finds climada

startup % climada_global exists afterwards
pwd % just to check where the job is running from
% N_pool_workers=24; % for parpool
% climada_global.parfor=1; % for parpool

% pool=parpool(N_pool_workers);

% define where data will be saved
output_folder=[climada_global.data_dir filesep 'isimip/results'];

clear params;
% parameters on cluster
params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
params.entity_prefix='FL1950';
params_damfun.filename_suffix='PAA1';
params_damfun.filepath=[climada_global.data_dir filesep 'isimip/entities/damfun'];
% local parameters (for testing only)
params.entity_folder=climada_global.entities_dir;
params.entity_prefix='';
params_damfun.filename_suffix='PAA1';
params_damfun.filepath='/Users/bguillod/Documents/work/ETH/floods/damage_functions/files_mine';

% generate a list of all countries to loop over
NatID_RegID=isimip_NatID_RegID; % get the mapping ISO3 - isimip country code
all_countries = NatID_RegID.ISO3;
% test: run only a subset of countries
all_countries = all_countries(1:4);

% do one file for each model combination
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
ghm = ghms{1};
forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};

for iGHM=1:length(ghms)
    for iForcing=1:length(forcings)
        ghm=ghms{iGHM};
        forcing = forcings{iForcing};
        
        for i=1:length(all_countries)
            country=all_countries(i);
            output_i = isimip2a_FL_countrydata_for_PIK(country, ghm, forcing, params, params_damfun);
            if i==1
                output_all = output_i;
            else
                output_all = cat(1, output_all, output_i);
            end
            
        end
        output_file = [output_folder filesep 'output_' ghm '_' forcing ''];
        csvwrite(output_file,output_all);
        
    end
end

% % climada_global.entities_dir='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   country='Switzerland';
%   ghm='CLM';
%   forcing='gswp3';
  
%   output_all2 = isimip2a_FL_countrydata_for_PIK(country, ghm, forcing, params, params_damfun);

% 
% params.grid_resolution='0150as';
% params.entity_prefix='FL1950_';
% isimip_gdp_entity('all',params,1950,2020);
%isimip_gdp_entity_TEST('all','0150as',1900,2018); % for parpool tests

%delete(pool)

exit % the cluster appreciates this, gives back memory etc.