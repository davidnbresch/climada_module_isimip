% batch job for cluster: bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_BGD_TS
% MODULE:
%   storm_europe
% NAME:
%   job_BGD_TS
% PURPOSE:
%   create Bangladesh probabilistic TS hazard set a caller
%   for climada_ts_hazard_set on the cluster
%
%   See PARAMETERS below before running
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/job_BGD_TS.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy data to cluster:      scp -r Documents/_GIT/climada_data/entities/BGD_Bangladesh_01x01.mat dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/entities/.
%   copy data to cluster:      scp -r Documents/_GIT/climada_data/tc_tracks/ibtracs/ibtracs.mat dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/tc_tracks/ibtracs/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_BGD_TS
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/_BGD_TS_ETOP.mat Documents/_GIT/climada_data/hazards/.
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/_BGD_TS_SRTM.mat Documents/_GIT/climada_data/hazards/.
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_BGD_TS
% EXAMPLE:
%   bsub -R "rusage[mem=500]" -n 24 matlab -nodisplay -singleCompThread -r job_BGD_TS
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20180101, copy from job_WISC
%-

% PARAMETERS
%
% one should only have to edit this section
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada
%scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
%%hazards_dir='/cluster/work/climate/dbresch/climada_data/hazards';


startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool


pool=parpool(N_pool_workers);

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
entity='BGD_Bangladesh_01x01.mat';
IBTrACS_ID='2007314N10093';
%interest_area=[]; % all centroids
interest_area=[89 92 21 23.5];
%
% the file with all IBTrACS alread read (see isimip_ibtracs_read, automatically invoked below)
ibtracs_save_file=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat'];


if exist(ibtracs_save_file,'file')
    load(ibtracs_save_file) % contains tc_track
else
    tc_track=isimip_ibtracs_read('all','',1,1);
end

% find the requested track
if ischar(IBTrACS_ID) % convert to integer for faster search
    IBTrACS_ID_no=str2double(strrep(strrep(IBTrACS_ID,'S','1'),'N','0'));
    fprintf('%s --> %i\n',IBTrACS_ID,IBTrACS_ID_no) % S->1, N->0
else
    IBTrACS_ID_no=IBTrACS_ID;
end

n_tracks=length(tc_track);IBTrACS_index=[];
for track_i=1:n_tracks
    if ~isempty(find(tc_track(track_i).ID_no==IBTrACS_ID_no, 1))
        IBTrACS_index=track_i;
    end
end % track_i

if ~isempty(IBTrACS_index)
    fprintf('processing %s (ID_no %i):\n',tc_track(IBTrACS_index).ID_str,IBTrACS_ID_no);
    
    % load the high-resolution entity
    entity=climada_entity_load(entity);
    if isempty(entity),return;end
    
    % restrict to area of interest (especially for time/memory resons)
    if ~isempty(interest_area)
        fprintf('> restricting to interest_area\n');
        interest_area_pos=(entity.assets.lon>interest_area(1) & entity.assets.lon<interest_area(2)) ...
            & (entity.assets.lat>interest_area(3) & entity.assets.lat<interest_area(4));
        entity.assets=climada_subarray(entity.assets,interest_area_pos);
    end
    
    % generate the TC windfield
    % to improve temporal resolution, consider climada_global.tc.default_min_TimeStep
    hazard_TC = climada_tc_hazard_set(tc_track(IBTrACS_index),'_BGD_Barisal_TC',entity);
    
    % generate the simply coarse-resolution (ETOPO)  TS surge field
    climada_ts_hazard_set(hazard_TC,'_BGD_Barisal_TS_ETOP');
    
    % generate the high-resolution (SRTM) TS surge field
    climada_ts_hazard_set(hazard_TC,'_BGD_Barisal_TS_SRTM','SRTM',0,1); % 0 for no check plot
    
else
    fprintf('ID_no %i not found in %s\n',IBTrACS_ID_no,ibtracs_save_file);
end

delete(pool)

exit % the cluster appreciates this, gives back memory etc.