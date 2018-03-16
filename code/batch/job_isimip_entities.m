% batch job for cluster: bsub -R "rusage[mem=5000]" -n 2 matlab -nodisplay -singleCompThread -r job_isimip_entities
% MODULE:
%   isimip
% NAME:
%   job_isimip_entities
% PURPOSE:
%   generate isimip entities for all countries
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
% David N. Bresch, dbresch@ethz.ch, 20180306, copy from job_isimip05
% David N. Bresch, dbresch@ethz.ch, 20180314, simple, not parallel
%-


% PARAMETERS
% one should only have to edit this section
%
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada

startup % climada_global exists afterwards
pwd % just to check where the job is running from
% N_pool_workers=24; % for parpool
% climada_global.parfor=1; % for parpool

% pool=parpool(N_pool_workers);

climada_global.entities_dir='/cluster/work/climate/dbresch/climada_data/isimip/entities';
params.grid_resolution='0150as';
params.entity_prefix='FL1950_';
isimip_gdp_entity('all',params,1950,2020);
%isimip_gdp_entity_TEST('all','0150as',1900,2018); % for parpool tests

%delete(pool)

exit % the cluster appreciates this, gives back memory etc.