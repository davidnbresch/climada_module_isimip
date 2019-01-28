% batch job for cluster: bsub -W 24:00 -R "rusage[mem=50000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip_calib_all
% MODULE:
%   isimip
% NAME:
%   job_isimip_calib_all
% PURPOSE:
%   calibrating shape and scale parameters for floods in isimip:
%   sample the parameter space (20 times 20) and do pattern search, with
%   different options
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r $climada_modules_folder/isimip/batch/jobs_isimip2a/job_isimip_calib_all.m bguillod@euler.ethz.ch:/cluster/home/bguillod/euler_jobs/
%   run on cluster:            bsub -R "rusage[mem=50000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip_calib_all
%       note: might need to request long wall time (e.g., -W 72:00) or to run in several batches by editing line 91 as follows:
%       regions_def=regions_def(1:15);
%       e.g. once replacing 1:15 with 1:4, then 5:8, etc.
%   copy results back local:   scp -r rsync -ar --include='*.csv' --include='*/' --exclude='*' bguillod@euler.ethz.ch:/cluster/work/climate/bguillod/climada_data/isimip/results/calibration $HOME/Documents/work/ETH/climada/DATA_OUTPUT/euler/isimip/results
%
% CALLING SEQUENCE:
%   bsub -W 24:00 -R "rusage[mem=50000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip_calib_all
% EXAMPLE:
%   bsub -W 24:00 -R "rusage[mem=50000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip_calib_all
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see code
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, copy from job_isimip_calib_01.m
%-


cd /cluster/home/bguillod/climada % to make sure the cluster finds climada
startup % climada_global exists afterwards
pwd % just to check where the job is running from

% all parameters to be defined
%(commented if the default value is the good one)


%% direct input
RegionID='EUR';
% years_range=[1992 2010];

%% params
params=struct;
% params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
% params.entity_prefix='FL1950';
% params.hazard_protection='flopros';
% params.entity_year=0;
% params.output_folder=[climada_global.data_dir filesep 'isimip/results/calibration'];
% params.keep_countries_0emdat=2;
% params.verbose_exclude_countries=0;

%% params_MDR
params_MDR=struct;
% params_MDR.remove_years_0YDS=0;
% params_MDR.pars_range={[0.0001,1],[0.0001,5]};
% params_MDR.damFun_xVals=0:0.5:15;

%% params_calibration
% either
params_calibration=struct;
params_calibration.calib_options.method='regular_sampling';
params_calibration.calib_options.params.n_per_dim=20;
% or
params_calibration.calib_options.method='patternsearch';
params_calibration.calib_options.params.random=0;
params_calibration.calib_options.params.nstart=5;
% params_calibration.calib_options.params.InitialMeshSize=0.25
% params_calibration.calib_options.params.step_tolerance=0.001;
% end or
params_calibration.type='dlog2';
params_calibration.underestimation_factor=2;
params_calibration.exclude_years_0totals=0; % also test 1 for dlog2
% params_calibration.MM_how='MMM';

%% params_computation
params_computation=struct;
% params_computation.years_range=[1971 2010];
% params_computation.do=0;





%% prepare to do a loop on RegionID and various parameters

% RegionID
fileID=fopen('/cluster/work/climate/bguillod/climada_data/isimip/RegID_names_isimip_flood.txt','r');
regions_def=strsplit(fscanf(fileID,'%c'),'\n');
fclose(fileID);
% skip PIS1 and PIS2 - not enough data there anyway
regions_def=regions_def(1:15);
regionIDs = {};
for i=1:length(regions_def)
temp = strsplit(regions_def{i},':');
regionsIDs{i} = temp{1};
end

% params_calibration, first regular_sampling, then patternsearch, both with
% dlog2 and RTarea, underestimation_factor=2
params_calibration_list={};
i=length(params_calibration_list)+1;
params_calibration_list{i}=struct;
params_calibration_list{i}.calib_options.method='regular_sampling';
params_calibration_list{i}.calib_options.params.n_per_dim=20;
params_calibration_list{i}.type='dlog2';
params_calibration_list{i}.underestimation_factor=2;
params_calibration_list{i}.exclude_years_0totals=0;
% now do dlog2 but removing years with 0 damages in total (either emdat or modelled)
i=length(params_calibration_list)+1;
params_calibration_list{i}=params_calibration_list{i-1};
params_calibration_list{i}.exclude_years_0totals=1;
% now do RTarea (not removing years with 0 damages in total as they're
% removed by definition)
i=length(params_calibration_list)+1;
params_calibration_list{i}=params_calibration_list{i-2};
params_calibration_list{i}.type='RTarea';

i=length(params_calibration_list)+1;
params_calibration_list{i}=struct;
params_calibration_list{i}.calib_options.method='patternsearch';
params_calibration_list{i}.calib_options.params.random=0;
params_calibration_list{i}.calib_options.params.nstart=3;
params_calibration_list{i}.type='dlog2';
params_calibration_list{i}.underestimation_factor=2;
params_calibration_list{i}.exclude_years_0totals=0;
% now do dlog2 but removing years with 0 damages in total (either emdat or modelled)
i=length(params_calibration_list)+1;
params_calibration_list{i}=params_calibration_list{i-1};
params_calibration_list{i}.exclude_years_0totals=1;
% now do RTarea
i=length(params_calibration_list)+1;
params_calibration_list{i}=params_calibration_list{i-2};
params_calibration_list{i}.type='RTarea';



%% function call for all
test=tic;
for i=1:length(regionsIDs)
    for j=1:length(params_calibration_list)
        fprintf('\n\n\n********** starting i=%s out of %s, j=%s out of %s **********\n\n\n',num2str(i),num2str(length(regionsIDs)),num2str(j),num2str(length(params_calibration_list)))
        [status,output_filename,output]=isimip_flood_calibration(regionsIDs{i},[1992 2010],params,params_MDR,params_calibration_list{j},params_computation);
    end
end
toc(test)

exit % the cluster appreciates this, gives back memory etc.
