% batch job for cluster: bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip03
% MODULE:
%   isimip
% NAME:
%   job_isimip03
% PURPOSE:
%   calibrating shape and scale parameters for floods in isimip.
%   THE REST OF THIS DOCUMENTATION HAS NOT BEEN UPDATED YET (TO DO).
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/euler_jobs/job_isimip03.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip03
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v03/.
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip03
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip03
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 2081030, copy from job_isimip03
%-

% PARAMETERS
% one should only have to edit this section
cd /cluster/home/bguillod/climada % to make sure the cluster finds climada
startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
parallel_what='optim'; % 'optim' to parallelize the optimization algorithm, 'parfor' to parallelize the computation of damages etc.


% prepare input parameters
RegionID='NAM';
years_range=[1990 2010];
%params
params=struct;
params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
params.hazard_protection = 'flopros';
params.subtract_matsiro = 0;
params.entity_year = 0;
%params_MDR
params_MDR=struct;
params_MDR.remove_years_0emdat=1;
params_MDR.remove_years_0YDS.do=0;
params_MDR.damFun_xVals=0:0.5:12;
%params_calibration
params_calibration=struct;
params_calibration.type='dlog2';
params_calibration.MM_how='MMM';
params_calibration.step_tolerance=0.05;

switch parallel_what
    case 'optim'
        climada_global.parfor=0;
        params_calibration.parallel=true;
    case 'parfor'
        climada_global.parfor=1;
        params_calibration.parallel=false;
    otherwise
        error('** unexpected value in parallel_what **');
end

pool=parpool(N_pool_workers);
[status,file_out]=isimip_flood_calibration(RegionID,years_range,params,params_MDR,params_calibration);
delete(pool)

if status
    fprintf('Calibration has succeeded\n');
    fprintf('Output saved in file %s\n',file_out);
else
    fprintf('Calibration has failed\n');
end


exit % the cluster appreciates this, gives back memory etc.
