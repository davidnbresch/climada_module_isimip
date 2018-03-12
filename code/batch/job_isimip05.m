% batch job for cluster: bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
% MODULE:
%   isimip
% NAME:
%   job_isimip05
% PURPOSE:
%   generate isimip probabilistic tropical cyclone (TC) hazard event sets based on Kerry
%   Emmanuel TC track files i.e. files such as Trial1_GB_dkmiroc_20thcal.mat
%
%   But use a coarse (1deg) grid instead of the (default 360as one)
%
%   At the end, copy the hazard event sets to dkrz (isimip)
%
%   Once all runs are done, consider using SPECIAL CODE2 at the end below
%   to create check plots
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/batch/job_isimip05.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/isimip/tc_tracks dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/isimip/.
%   copy single data to cluster:scp -r Documents/_GIT/climada_data/isimip/tc_tracks/Trial3_GB_dkgfdl_piControlcal dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_data/isimip/tc_tracks/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v04/.
%   copy results to dkrz:      scp -r /cluster/work/climate/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
%
%   run_on_desktop=1; % to test the job on a desktop
%   job_isimip05
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling job_isimip05 (all in
%       MATLAB command window), it sets the number of parallel pool workers 
%       to two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.  
%       Default=0 for cluster.
% OUTPUTS:
%   to disk, see PARAMETERS and climada folders
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20171123, copy from job_isimip04
% David N. Bresch, dbresch@ethz.ch, 20180309, switched to /cluster/work/climate/dbresch
% David N. Bresch, dbresch@ethz.ch, 20180310, desktop option added (run_on_desktop)
%-

% PARAMETERS
%
FAST_TEST=0; % default=0, if =1, set -R "rusage[mem=500]"
%
cluster_climada_root_dir='/cluster/home/dbresch/climada'; % to make sure the cluster finds climada
cluster_N_pool_workers=24; % number of parpool workers on pool (same as argument in bsub -n ..)
desktop_N_pool_workers= 2; % number of parpool workers on desktop
%
% the list of TC track files to be processed (see SPECIAL CODE below)
track_files={
    'Trial1_GB_dkmiroc_20thcal'
    'Trial1_GB_dkmiroc_rcp60cal'
    };


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
climada_global.parfor=1; % for parpool
t0=clock;
% eee: end of admin (do not edit until here)


climada_global.tc.default_min_TimeStep=0.25; % 15 min

% prepare centroids (full globe leads to segmentation fault)
centroids=climada_centroids_load('GLB_NatID_grid_0360as_adv_1');
 
% sub-select only the regular grid points
% ---------------------------------------
%centroids.lon=round(centroids.lon,4);centroids.lat=round(centroids.lat,4); % as we have a precision issue
%pos= centroids.lon==ceil(centroids.lon) & centroids.lat==ceil(centroids.lat);
grid_pos=find(centroids.centroid_ID>3e6);
% note that we added +0.05 when we generated the grid to set it apart
centroids.lon              =centroids.lon(grid_pos)-0.05; % subtract bias
centroids.lat              =centroids.lat(grid_pos)-0.05; % subtract bias
centroids.centroid_ID      =centroids.centroid_ID(grid_pos);
centroids.distance2coast_km=centroids.distance2coast_km(grid_pos);
if isfield(centroids,'NatID'),centroids=rmfield(centroids,'NatID');end
if isfield(centroids,'admin1_ID'),centroids=rmfield(centroids,'admin1_ID');end
if isfield(centroids,'ISO3_list'),centroids=rmfield(centroids,'ISO3_list');end

if FAST_TEST
    fprintf('\n !!! FAST TEST mode - only a small subset of centroids and tracks !!!\n\n');
    % reduce centroids for TEST
    centroids.lon=centroids.lon(1:1000:end);
    centroids.lat=centroids.lat(1:1000:end);
    centroids.centroid_ID=centroids.centroid_ID(1:1000:end);
    centroids.distance2coast_km=centroids.distance2coast_km(1:1000:end);
    centroids.NatID=centroids.centroid_ID;
end

%for file_i=1:length(track_files)
for file_i=2:length(track_files)
    
    hazard_name        = [track_files{file_i} '_1deg_prob'];
    hazard_set_file    = [climada_global.hazards_dir filesep hazard_name];
    tc_track_prob_name = [climada_global.data_dir    filesep 'isimip' filesep 'tc_tracks' filesep track_files{file_i} '_tracks_prob.mat'];
    
    if exist(tc_track_prob_name,'file')
        fprintf('loading from %s\n',tc_track_prob_name);
        load(tc_track_prob_name) % 20180309
    else
        tc_track=isimip_tc_track_load(track_files{file_i},'both',180,-1); % both hemispheres
        
        % create probabilistic tracks
        [~,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,0); % wind speed decay at track nodes after landfall
        tc_track  = climada_tc_random_walk(tc_track); % overwrites tc_track to save memory
        tc_track  = climada_tc_track_wind_decay(tc_track, p_rel,0); % add the inland decay correction to all probabilistic nodes
        
        fprintf('storing prob tracks %s\n',tc_track_prob_name)
        save(tc_track_prob_name,'tc_track',climada_global.save_file_version);
    end
    
    %     if check_plots
    %         % plot the tracks
    %         figure('Name','TC tracks','Color',[1 1 1]); hold on
    %         for event_i=1:length(tc_track) % plot all tracks
    %             plot(tc_track(event_i).lon,tc_track(event_i).lat,'-b');
    %         end % event_i
    %         % overlay historic (to make them visible, too)
    %         for event_i=1:length(tc_track)
    %             if tc_track(event_i).orig_event_flag
    %                 plot(tc_track(event_i).lon,tc_track(event_i).lat,'-r');
    %             end
    %         end % event_i
    %     end % check_plots
    
    n_track_tot=length(tc_track);
    n_track_end_1=floor(n_track_tot/2);
    fprintf('running in two lumps: tracks 1..%i, %i..%i\n',n_track_end_1,n_track_end_1+1,n_track_tot);
    
    %for ii=1:2
   %     
    %    if ii==1
    %        tc_track=tc_track(1:n_track_end_1);
    %        hazard_set_file    = [climada_global.hazards_dir filesep hazard_name '_1'];
    %    elseif ii==2
    %        clear tc_track
    %        fprintf('loading from %s\n',tc_track_prob_name);
    %        load(tc_track_prob_name)
            tc_track=tc_track(n_track_end_1+1:end);
            hazard_set_file    = [climada_global.hazards_dir filesep hazard_name '_2'];
     %   end
        
        if FAST_TEST,tc_track=tc_track(1:1000);end % small subset for TEST
        fprintf('using %i centroids, %i (prob) tracks\n',length(centroids.lon),length(tc_track))
        
        isimip_tc_hazard_set(tc_track,hazard_set_file,centroids,0,hazard_name);
   % end % ii
   
   % if you slit in two, use hazard=climada_hazard_merge(hazard,hazard_2,'events');
    
end % file_i

% copy results to dkrz (no closing ; to log success):
% ---------------------
%[status,result]=system('scp -r    /cluster/work/climate/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')
%[status,result] =system('scp -r -v /cluster/work/climate/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')

% %
% % SPECIAL CODE2 to inspect results
% % ----------------------------------
% result_hazard_dir='/Users/bresch/polybox/isimip/';
% fig_dir='/Users/bresch/Desktop/isimip';fig_ext='png';
% if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
% params.figure_scale=0; % no geographical scale on figure
% params.blue_ocean=1;   % no geographical scale on figure
% n_files=length(track_files);
% for file_i=1:n_files % define them as above in PARAMETERS
%     fprintf('processing %s: (%i of %i)\n',track_files{file_i},file_i,n_files)
%     hazard=climada_hazard_load([result_hazard_dir filesep track_files{file_i}]);
%     figure;climada_hazard_plot(hazard,0,[],params); % plot max intensity
%     fN=strrep(track_files{file_i},'.mat','');xlim([-180 180]);
%     ylim1=ylim;if ylim1(2)>90,ylim1(2)=90;end;if ylim1(2)<-90,ylim1(2)=-90;end
%     title_str=strrep(fN,'_',' ');ylim(ylim1);
%     title([title_str ' (max intensity)']);
%     saveas(gcf,[fig_dir filesep fN '.' fig_ext],fig_ext);
%     close all
% end % file_i


if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit