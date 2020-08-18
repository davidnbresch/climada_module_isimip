% batch job for cluster: bsub -W 48:00 -R "rusage[mem=9000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip06a
% MODULE:
%   isimip
% NAME:
%   job_isimip06a
% PURPOSE:
%   generate isimip tropical cyclone (TC) hazard event sets based on Kerry
%   Emmanuel TC track files as sent by Tobias Geiger 24.6.2020, see track_files={ ... below
%
%   Please consider running SPECIAL CODE to sort track files by size
%   manually (see below). At the end, the batch job produces check figures
%   (see SPECIAL CODE2) and tries to copy the hazard event sets to DKRZ
%   (isimip).
%
%   Note: Currently, June 2020, the automatic copy to DKRZ does not work for
%   security reasons, hence issue the rcp command manually, code commented
%   out below.
%
%   See PARAMETERS before running this batch job
%
%   NOTE: on Euler, we moved all previous /cluster/work/climate/dbresch/climada_data/hazards/Trial*.mat
%         hazards to ../hazards/Trials_old in order to avoid overwrite etc.
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/batch/job_isimip06a.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy single tc_track data to cluster: scp -r Documents/_GIT/climada_data/isimip/tc_tracks/* dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/isimip/tc_tracks/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip06a
%   check progress:            ls -la /cluster/work/climate/dbresch/climada_data/hazards/Trial*.mat
%   check with bjobs and at the end check the file /cluster/home/dbresch/euler_jobs/lsf.* (stdout)
%   generall, see http://www.clusterwiki.ethz.ch/brutus/LSF_mini_reference
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/Trial*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v04/.
%   copy results to dkrz:      scp -r /cluster/work/climate/dbresch/climada_data/hazards/Trial*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.
%
%   other option, a LSF pool, see http://www.clusterwiki.ethz.ch/brutus/Parallel_MATLAB_and_Brutus
%
% JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
% 128646395  dbresch RUN   bigmem.120 eu-login-10 24*eu-g1-03 *isimip06a Jun 25 11:00 % see below
%
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip06a
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip06a
%
%   run_on_desktop=1; % to test the job on a desktop
%   job_isimip06a
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling job_isimip06a (all in
%       MATLAB command window), it sets the number of parallel pool workers
%       to two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.
%       Default=0 for cluster.
% OUTPUTS:
%   to disk, see PARAMETERS and climada folders
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20200625, copy from job_isimip04t4
% David N. Bresch, dbresch@ethz.ch, 20200625, run with FAST_TEST=1 all fine
%-


% PARAMETERS
%
FAST_TEST=0; % default=0, if =1, set -R "rusage[mem=500]"
%
cluster_climada_root_dir='/cluster/home/dbresch/climada'; % to make sure the cluster finds climada
cluster_N_pool_workers=24; % number of parpool workers on pool (same as argument in bsub -n ..)
desktop_N_pool_workers= 2; % number of parpool workers on desktop
%
% the list of TC track files to be processed (see SPECIAL CODE below, sorted in ascending size)
track_files={
    'Trial1_GB_dkipsl_20thcal' % started Jun 25 11:00
    'Trial1_GB_dkmiroc_rcp85cal' % done Jun 25 20:15
    'Trial1_GB_dkgfdl_rcp85cal'  % done Jun 26 01:53
    'Trial1_GB_dkhad_rcp85cal'   % done Jun 26 06:50
    'Trial1_GB_dkipsl_rcp85cal'  % done Jun 26 12:56
    'Trial1_GB_dkipsl_rcp60cal'  % done Jun 26 19:16
    'Trial1_GB_dkipsl_rcp26cal'  % done  Jun 27 01:30
    'Trial4_GB_dkipsl_20thcal'   % done Jun 27 07:37, job then failed (segmentation fault)
    'Trial4_GB_dkipsl_rcp26cal'   % done Jun 28 06:22 
    };
%
% % SPECIAL CODE to sort track files by size (run this on command line to
% % obtain above list). Ascending size helps to run until we encounter
% % memory problems and can then tackle them (if needed)
% % ---------------------------------------------------------------------
% dd_name={};dd_bytes=[]; % init
% tc_tracks_dir=[climada_global.data_dir filesep 'isimip/tc_tracks'];
% dd=dir(tc_tracks_dir);
% for i=1:length(dd)
%     if dd(i).isdir && length(dd(i).name)>2
%         check_file=[tc_tracks_dir filesep dd(i).name filesep dd(i).name '.mat'];
%         ddd=dir(check_file);
%         dd_name{end+1}=dd(i).name;
%         dd_bytes(end+1)=ddd(1).bytes;
%     end % dd(i).isdir
% end % i
% [~,pos] = sort(dd_bytes);
% dd_bytes=dd_bytes(pos);
% dd_name=dd_name(pos);
% fprintf('track files sorted ascening by size (%2.2g..%2.2g bytes):\n',dd_bytes(1),dd_bytes(end));
% for i=1:length(dd_name),fprintf('''%s''\n',dd_name{i});end
% fprintf('--> copy paste this into track_files in %s\n','job_isimip06a');


% AAA: some admin to start with (up to % EEE standard code)
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
if exist('startupp','file'),startupp;end % sets data folders on cluster
fprintf('executing in %s\n',pwd) % just to check where the job is running from
climada_global.parfor=1; % for parpool
t0=clock;
% EEE: end of admin (do not edit until here)


climada_global.tc.default_min_TimeStep=0.25; % in hours, i.e. 15 min

% prepare centroids (full globe leads to segmentation fault)
centroids_S=climada_centroids_load('GLB_NatID_grid_0360as_adv_1');

if FAST_TEST
    fprintf('\n !!! FAST TEST mode - only a small subset of centroids and tracks !!!\n\n');
    % reduce centroids for TEST
    centroids_S.lon=centroids_S.lon(1:1000:end);
    centroids_S.lat=centroids_S.lat(1:1000:end);
    centroids_S.centroid_ID=centroids_S.centroid_ID(1:1000:end);
    centroids_S.distance2coast_km=centroids_S.distance2coast_km(1:1000:end);
    centroids_S.NatID=centroids_S.centroid_ID;
end

centroids_N = centroids_S; % Northern hemisphere
lat_pos=find(centroids_N.lat>0 & centroids_N.lat<=60);
centroids_N.lon=centroids_N.lon(lat_pos);
centroids_N.lat=centroids_N.lat(lat_pos);
centroids_N.centroid_ID=centroids_N.centroid_ID(lat_pos);
centroids_N.distance2coast_km=centroids_N.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_N.NatID))=[];
centroids_N.NatID=centroids_N.NatID(lat_pos);

% centroids_NE = centroids_N; % Northeastern hemisphere
% lon_pos=find(centroids_NE.lon>=0);
% centroids_NE.lon=centroids_NE.lon(lon_pos);
% centroids_NE.lat=centroids_NE.lat(lon_pos);
% centroids_NE.centroid_ID=centroids_NE.centroid_ID(lon_pos);
% centroids_NE.distance2coast_km=centroids_NE.distance2coast_km(lon_pos);
% lon_pos(lon_pos>length(centroids_NE.NatID))=[];
% centroids_NE.NatID=centroids_NE.NatID(lon_pos);
%
% centroids_NW = centroids_N; % Northwestern hemisphere
% lon_pos=find(centroids_NW.lon<0);
% centroids_NW.lon=centroids_NW.lon(lon_pos);
% centroids_NW.lat=centroids_NW.lat(lon_pos);
% centroids_NW.centroid_ID=centroids_NW.centroid_ID(lon_pos);
% centroids_NW.distance2coast_km=centroids_NW.distance2coast_km(lon_pos);
% lon_pos(lon_pos>length(centroids_NW.NatID))=[];
% centroids_NW.NatID=centroids_NW.NatID(lon_pos);

lat_pos=find(centroids_S.lat<=0 & centroids_S.lat>=-60); % Southern hemisphere
centroids_S.lon=centroids_S.lon(lat_pos);
centroids_S.lat=centroids_S.lat(lat_pos);
centroids_S.centroid_ID=centroids_S.centroid_ID(lat_pos);
centroids_S.distance2coast_km=centroids_S.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_S.NatID))=[];
centroids_S.NatID=centroids_S.NatID(lat_pos);

% centroids_SE = centroids_S; % Southeastern hemisphere
% lon_pos=find(centroids_SE.lon>=0);
% centroids_SE.lon=centroids_SE.lon(lon_pos);
% centroids_SE.lat=centroids_SE.lat(lon_pos);
% centroids_SE.centroid_ID=centroids_SE.centroid_ID(lon_pos);
% centroids_SE.distance2coast_km=centroids_SE.distance2coast_km(lon_pos);
% lon_pos(lon_pos>length(centroids_SE.NatID))=[];
% centroids_SE.NatID=centroids_SE.NatID(lon_pos);
%
% centroids_SW = centroids_S; % Southwestern hemisphere
% lon_pos=find(centroids_SW.lon<0);
% centroids_SW.lon=centroids_SW.lon(lon_pos);
% centroids_SW.lat=centroids_SW.lat(lon_pos);
% centroids_SW.centroid_ID=centroids_SW.centroid_ID(lon_pos);
% centroids_SW.distance2coast_km=centroids_SW.distance2coast_km(lon_pos);
% lon_pos(lon_pos>length(centroids_SW.NatID))=[];
% centroids_SW.NatID=centroids_SW.NatID(lon_pos);

for file_i=1:length(track_files)
    
    tc_track=isimip_tc_track_load(track_files{file_i},'N',180,-1); % Northern hemisphere
    if FAST_TEST,tc_track=tc_track(1:100);end % small subset for TEST
    hazard_name=[track_files{file_i} '_N_0360as'];
    %hazard_set_file=[scratch_dir filesep hazard_name];
    %isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_N,0,hazard_name);
    isimip_tc_hazard_set(tc_track,hazard_name,centroids_N,0,hazard_name);
    
    tc_track=isimip_tc_track_load(track_files{file_i},'S',180,-1); % Southern hemisphere
    if FAST_TEST,tc_track=tc_track(1:100);end % small subset for TEST
    hazard_name=[track_files{file_i} '_S_0360as'];
    %hazard_set_file=[scratch_dir filesep hazard_name];
    %isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_S,0,hazard_name);
    isimip_tc_hazard_set(tc_track,hazard_name,centroids_S,0,hazard_name);
    
    clear tc_track % might help with memory usage
    
end % file_i


% copy results to dkrz (no closing ; to log success): commented 20200625,
% as no keyaccess allowed on DKRZ temporarily for security reasons
% ---------------------
%[status,result]=system('scp -r -v /cluster/work/climate/dbresch/climada_data/hazards/Trial*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')


% % SPECIAL CODE2 to inspect results
% % --------------------------------
%define track_files as in PARAMETERS above
result_hazard_dir=climada_global.hazards_dir;
fig_dir=[climada_global.results_dir filesep 'isimip'];fig_ext='png';
%result_hazard_dir='/Users/bresch/polybox/isimip/hazards';
%fig_dir='/Users/bresch/Desktop/isimip';fig_ext='png';
if ~isfolder(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
params.figure_scale=0; % no geographical scale on figure
params.blue_ocean=1; % no geographical scale on figure
marker_size=2;
figure_visible='off'; % ='off' or 'on', 'off' on cluster
n_track_files=length(track_files);
for file_i=1:n_track_files
    fprintf('processing %s_{N|S}_0360as: (%i of %i)\n',track_files{file_i},file_i,n_track_files)
    for NS_i=1:2
        if NS_i==1
            hazard_name=[track_files{file_i} '_N_0360as'];
        else
            hazard_name=[track_files{file_i} '_S_0360as'];
        end
        hazard=climada_hazard_load([result_hazard_dir filesep hazard_name]);
        plt_fig=figure('Name','isimip06a','NumberTitle','off','Position',[9 824 960 520],'Visible',figure_visible);
        res=climada_hazard_plot(hazard,0,marker_size,params); % plot max intensity
        % res.X and .Y hold coord, .VALUE holds intensity
        xlim([-180 180]);
        ylim1=ylim;if ylim1(2)>90,ylim1(2)=90;end;if ylim1(2)<-90,ylim1(2)=-90;end;ylim(ylim1);
        title_str=strrep(hazard_name,'_',' ');title([title_str ' (max intensity)']);
        saveas(plt_fig,[fig_dir filesep hazard_name '.' fig_ext],fig_ext);
        close all % to clear memory, mainly
    end % NS_i
end % file_i
fprintf('figures written to %s\n',fig_dir)

if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit
