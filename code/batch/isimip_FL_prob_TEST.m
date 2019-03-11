%function res=isimip_FL_prob_TEST
% climada isimip_FL_prob_TEST
% MODULE:
%   _LOCAL
% NAME:
%   isimip_FL_prob_TEST
% PURPOSE:
%   batch job to TEST flood hazard in one country
%
%   see PARAMETERS
%   for speedup, consider climada_global.parfor=1
%   to run this on the cluster as batch, see batch_job_template
%   last time run on cluster on command line (not recommended)
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%
%   copy job to cluster: scp -r Documents/_GIT/climada_modules/_LOCAL/code/isimip_FL_prob_TEST.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/. isimip_FL_prob_TEST.m         
%   copy data from cluster: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/GBR_FL_test* /Users/bresch/Documents/_GIT/climada_data/results/GBR_FL_test/.
%
%   for more, see http://www.clusterwiki.ethz.ch/brutus/Parallel_MATLAB_and_Brutus
% CALLING SEQUENCE:
%   isimip_FL_prob_TEST
% EXAMPLE:
%   bsub -W 2:00 -R "rusage[mem=1000]" -n 2 matlab -nodisplay -singleCompThread -r GBR_FL_test
%   -W: time, here 2 hours
%   mem: memory, for large jobs, request e.g. 9000
%   -n: number of cluster workers, here 2
%
%   run_on_desktop=1; % to test the job on a desktop
%   isimip_FL_prob_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling GBR_FL_test (all in
%       MATLAB command window), it sets the number of parallel pool workers
%       to two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.
%       Default=0 for cluster.
%   GBR_FL_test
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20190201, initial
% David N. Bresch, david.bresch@gmail.com, 20190202, fit for cluster
% David N. Bresch, david.bresch@gmail.com, 20190311, from CHN_FL_test to isimip_FL_prob_TEST
%-


% CLUSTER JOB PARAMETERS and setup - see PARAMETERS below for specifics
%
FAST_TEST=0; % default=0, if =1, set -R "rusage[mem=500]"
%
cluster_climada_root_dir='/cluster/home/dbresch/climada'; % to make sure the cluster finds climada
cluster_N_pool_workers= 2; % number of parpool workers on pool (same as argument in bsub -n ..)
desktop_N_pool_workers= 2; % number of parpool workers on desktop

% aaa: some admin to start with (up to % eee standard code)
if ~exist('run_on_desktop','var'),run_on_desktop=[];end
if isempty(run_on_desktop),run_on_desktop=0;end % default=0, =1 to run job on mac
if run_on_desktop % for parpool on desktop
    N_pool_workers=desktop_N_pool_workers;
    pool_active=gcp('nocreate');
    if isempty(pool_active),pool=parpool(N_pool_workers);end
else
    cd(cluster_climada_root_dir)
    N_pool_workers=cluster_N_pool_workers;
    pool=parpool(N_pool_workers);
end
startup % climada_global exists afterwards
if exist('startupp','file'),startupp;end
fprintf('executing in %s\n',pwd) % just to check where the job is running from
climada_global.parfor=1; % for parpool, see e.g. climada_tc_hazard_set
t0=clock;
% eee: end of admin (do not edit until here)


% PARAMETERS
%
ISO3='GBR'; % the country ISO3-code, such as 'GBR' for Great Britain
%
% the folder with isimip data on the cluster
cluster_data_folder='/cluster/work/climate/dbresch/climada_data/isimip/';
%
models={
    'FLOOD_CLM_gswp3'
    'FLOOD_CLM_princeton'
    'FLOOD_CLM_watch'
    'FLOOD_CLM_wfdei'
    'FLOOD_DBH_gswp3'
    'FLOOD_DBH_princeton'
    'FLOOD_DBH_watch'
    'FLOOD_DBH_wfdei'
    'FLOOD_H08_gswp3'
    'FLOOD_H08_princeton'
    'FLOOD_H08_watch'
    'FLOOD_H08_wfdei'
    'FLOOD_JULES-TUC_gswp3'
    'FLOOD_JULES-TUC_princeton'
    'FLOOD_JULES-TUC_wfdei'
    'FLOOD_JULES-UoE_gswp3'
    'FLOOD_JULES-UoE_princeton'
    'FLOOD_JULES-UoE_watch'
    'FLOOD_LPJmL_gswp3'
    'FLOOD_LPJmL_princeton'
    'FLOOD_LPJmL_watch'
    'FLOOD_LPJmL_wfdei'
    'FLOOD_MATSIRO_gswp3'
    'FLOOD_MATSIRO_princeton'
    'FLOOD_MATSIRO_watch'
    'FLOOD_MPI-HM_gswp3'
    'FLOOD_MPI-HM_princeton'
    'FLOOD_MPI-HM_watch'
    'FLOOD_MPI-HM_wfdei'
    'FLOOD_ORCHIDEE_gswp3'
    'FLOOD_ORCHIDEE_princeton'
    'FLOOD_ORCHIDEE_watch'
    'FLOOD_ORCHIDEE_wfdei'
    'FLOOD_PCR-GLOBWB_gswp3'
    'FLOOD_PCR-GLOBWB_princeton'
    'FLOOD_PCR-GLOBWB_watch'
    'FLOOD_PCR-GLOBWB_wfdei'
    'FLOOD_VEGAS_gswp3'
    'FLOOD_VIC_gswp3'
    'FLOOD_VIC_princeton'
    'FLOOD_VIC_watch'
    'FLOOD_VIC_wfdei'
    'FLOOD_WaterGAP_gswp3'
    'FLOOD_WaterGAP_princeton'
    'FLOOD_WaterGAP_watch'
    'FLOOD_WaterGAP_wfdei'
    };
%
if FAST_TEST,models=models(1:3);end % only first three models
%
ext_nomg=[      '_0_gev_0.1_FL1950_' ISO3 '_0150as_mFRCmatsiro_FL.mat'];
ext_mngd=['_flopros_gev_0.1_FL1950_' ISO3 '_0150as_mFRCmatsiro_FL.mat'];
%
% explicit list of all models (just FYI):
% FLOOD_CLM_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_CLM_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_DBH_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_H08_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-TUC_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_JULES-UoE_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_LPJmL_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MATSIRO_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_MPI-HM_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_ORCHIDEE_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_PCR-GLOBWB_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VEGAS_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VEGAS_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VEGAS_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_VIC_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_gswp3_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_gswp3_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_gswp3_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_princeton_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_princeton_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_princeton_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_watch_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_watch_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_watch_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_wfdei_0_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_wfdei_100_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat
% FLOOD_WaterGAP_wfdei_flopros_gev_0.1_FL1950_GBR_0150as_mFRCmatsiro_FL.mat


% load entity and add default damage functions (inlcudes one for FL)
entity=climada_entity_load([cluster_data_folder filesep 'entities' filesep 'FL1950_' ISO3 '_0150as_entity_dnb'],1); % no-save
entity_temp=climada_entity_load('entity_template'); % load template damage functions
entity.damagefunctions=entity_temp.damagefunctions; % replace to include FL damage  function

% % now adjust hand-made
% FL_dmf=climada_damagefunctions_plot(entity,'FL',1); % get default FL curve
% FL_dmf.PAA=[0 .01 .02 .05 .10 .20 .30 .40 .50 .60 .65 .7 .75 .8 .85 .875 .9 .925 .95];
% FL_dmf.MDD=[0 .01 .02 .05 .10 .20 .30 .40 .50 .60 .65 .675 .70 .71 .72 .73 .74 .75 .76];
% %FL_dmf.MDD=[0 .001  .01  .02  .05  .1   .20  .30 .40 .50 .6 .65 .70 .71 .72 .73 .74 .75 .76];
% %FL_dmf.MDD=[0 .0001 .001 .005 .01  .05  .1   .2  .3  .4  .5 .6  .65 .7  .75 .8  .85 .9  .9];
% %FL_dmf.MDD=[0 .0001 .001 .002 .003 .004 .005 .1  .2  .3  .4 .45 .5  .5  .5  .5  .5  .5  .5];
% FL_dmf.MDR=FL_dmf.PAA.*FL_dmf.MDD;
% entity=climada_damagefunctions_replace(entity,FL_dmf); % replace with new

% init output
n_models=length(models);
res.model_i=[];
res.event_ID=[];
res.damage_nomg=[];
res.damage_mngd=[];
res.frequency_nomg=[];
res.frequency_mngd=[];
res.ED=zeros(2,n_models);

for model_i=1:n_models
    
    fprintf('- dealing with model %s (%i of %i)\n',models{model_i},model_i,n_models);
    
    % load hazard with no flood protection (kind of no management, nomg)
    hazard_nomg=climada_hazard_load([cluster_data_folder filesep 'hazards' filesep models{model_i} ext_nomg]);
    
    if ~isfield(entity.assets,'hazard') % encode entity first time (takes time)
        entity=climada_assets_encode(entity,hazard_nomg);
        fprintf('saving encoded entity as %s ..',entity.assets.filename);
        save(entity.assets.filename,'entity')
        fprintf(' done\n');
    end
    
    EDS=climada_EDS_calc(entity,hazard_nomg);
    EDS(1).annotation_name='flood without management';
    
    % load hazard with flood protection (kind of managed, mngd)
    hazard_mngd=climada_hazard_load([cluster_data_folder filesep 'hazards' filesep models{model_i} ext_mngd]);
    
    %hazard_mngd=hazard; % generate hazard with mock FL management
    %pos=find(hazard_mngd.intensity>0);
    %hazard_mngd.intensity(pos)=max(hazard_mngd.intensity(pos)-1,0); % subtract 1m flood height
    EDS(2)=climada_EDS_calc(entity,hazard_mngd);
    EDS(2).annotation_name='flood with management';
    fprintf('  approx effect of management: %f%%\n',(EDS(2).ED-EDS(1).ED)/EDS(1).ED*100)
    % approx effect of management: -68.923298%
    % with FLOOD_VIC_princeton_100... instead: approx effect of management: 356.969439%
    
    fig1=figure('Visible','off');climada_EDS_DFC(EDS);
    saveas(fig1,[climada_global.results_dir filesep ISO3 '_FL_test_DFC_' models{model_i}],'png'); % in /cluster/work/climate/dbresch/climada_data/results
    delete(fig1)
    
    % append results
    res.model_i =[res.model_i  ones(1,hazard_nomg.event_count)*model_i];
    res.event_ID=[res.event_ID hazard_nomg.event_ID];
    res.damage_nomg=[res.damage_nomg EDS(1).damage];
    res.damage_mngd=[res.damage_mngd EDS(2).damage];
    res.frequency_nomg=[res.frequency_nomg EDS(1).frequency];
    res.frequency_mngd=[res.frequency_mngd EDS(2).frequency];
    
    res.ED(1,model_i)=EDS(1).ED;
    res.ED(2,model_i)=EDS(2).ED;
    
end % model_i

res.models=models;
res.n_years=length(res.damage_nomg); % events are just years

fprintf('dealt with %i models and %i years in total\n',n_models,n_years);
res.EDS = EDS; % to fill fields
res.EDS(1).damage=res.damage_nomg;
res.EDS(1).frequency=(1:res.n_years)*0+1/res.n_years;
res.EDS(1).ED_at_centroid=[];
res.EDS(1).event_ID=[];
res.EDS(1).ED=res.EDS(1).damage*res.EDS(1).frequency';
res.EDS(2).damage=res.damage_mngd;
res.EDS(2).frequency=res.EDS(1).frequency;
res.EDS(2).ED_at_centroid=[];
res.EDS(2).event_ID=[];
res.EDS(2).ED=res.EDS(2).damage*res.EDS(2).frequency';

fig1=figure('Visible','off');climada_EDS_DFC(res.EDS);
saveas(fig1,[climada_global.results_dir filesep ISO3 '_FL_test_DFC_combined'],'png'); % in /cluster/work/climate/dbresch/climada_data/results
delete(fig1)

save_filename=[climada_global.results_dir filesep ISO3 '_FL_test_data'];
fprintf('saving result data as %s\n',save_filename);
save(save_filename,'res');

% get EM-DAT for ISO3 FL (and scale to USD bn, as standard in isimip)
em_data=emdat_read('',ISO3,'FL',1,1);
em_data.damage          = em_data.damage          /entity.assets.currency_unit;
em_data.damage_orig     = em_data.damage_orig     /entity.assets.currency_unit;
em_data.DFC.damage      = em_data.DFC.damage      /entity.assets.currency_unit;
em_data.DFC_orig.damage = em_data.DFC_orig.damage /entity.assets.currency_unit;

fig1=figure('Visible','off');
[~,~,legend_str,legend_handle]=climada_EDS_DFC(res.EDS);
[legend_str,legend_handle]=emdat_barplot(em_data,'','','EM-DAT indexed',legend_str,legend_handle,'southeast');
xlim([0 250]);
saveas(fig1,[climada_global.results_dir filesep ISO3 '_FL_test_DFC_combined_emdat'],'png'); % in /cluster/work/climate/dbresch/climada_data/results
delete(fig1)

% special code to construct one hazard set with all simulations:
% --------------------------------------------------------------
fprintf('\nconstruct one hazard set with all %i simulations:\n',n_models);
n_models=length(models);
info.n_models=n_models;
info.models=models;
info.ext_mngd=ext_mngd;
info.model_i=[];
info.event_ID=[];
info.frequency=[];
hazard=[];
for model_i=1:n_models
    fprintf('- dealing with model %s (%i of %i)\n',models{model_i},model_i,n_models);
    hazard_mngd=climada_hazard_load([cluster_data_folder filesep 'hazards' filesep models{model_i} ext_mngd]);
    info.model_i = [info.model_i   ones(1,hazard_mngd.event_count)*model_i];
    info.event_ID= [info.event_ID  hazard_mngd.event_ID];
    info.frequency=[info.frequency hazard_mngd.frequency];
    hazard_mngd.yyyy   =str2num(hazard_mngd.yyyy)';
    hazard_mngd.mm     =str2num(hazard_mngd.mm)';
    hazard_mngd.dd     =str2num(hazard_mngd.dd)';
    hazard_mngd.datenum=hazard_mngd.datenum';
    if isempty(hazard)
        hazard=hazard_mngd;
    else
        hazard=climada_hazard_merge(hazard,hazard_mngd,'events');
    end
end % model_i
% re-define frequency
hazard.orig_years=hazard.orig_event_count; % we know each event is a year
hazard.frequency=hazard.frequency*0+1/hazard.orig_years;
hazard.comment='a collection of isimip FL hazars sets, see hazard.info';
fprintf('NOTE: event frequency redefined as 1/%i years (%i models combined)\n',hazard.orig_years,n_models);
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity); % update
hazard.info=info; % add info
hazard_filename=[climada_global.results_dir filesep ISO3 '_FL_test_merged'];
fprintf('saving combined hazard as %s ..',hazard_filename);
save(hazard_filename,'hazard');
fprintf(' done\n');

fprintf('job execution took %f sec.\n',etime(clock, t0))
if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit

% % and to assess a ptf locally (after copying {ISO3}_FL_test_merged.mat to local)
% entity=climada_entity_load('GBR_UnitedKingdom_10x10');
% hazard=climada_hazard_load([climada_global.results_dir  filesep 'GBR_FL_test' filesep 'GBR_FL_test_merged']);
% entity=climada_assets_encode(entity,hazard);
% EDS=climada_EDS_calc(entity,hazard);