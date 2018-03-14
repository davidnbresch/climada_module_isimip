% batch job for cluster: bsub -W 8:00 -R "rusage[mem=9000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip04
% MODULE:
%   isimip
% NAME:
%   job_isimip04
% PURPOSE:
%   generate isimip tropical cyclone (TC) hazard event sets based on Kerry
%   Emmanuel TC track files for 4th (corrected) batch, i.e. files such as
%   Trial4_GB_dkgfdl_20thcal.mat etc.
%
%   At the end, copy the hazard event sets to dkrz (isimip)
%
%   See PARAMETERS before running this and run SPECIAL CODE (also below) in
%   order to prepare the lostz of track files in ascending size order
%
%   Once all runs are done, consider using SPECIAL CODE2 at the end below
%   to create check plots
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/batch/job_isimip04.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   check progress:            ls -la /cluster/work/climate/dbresch/climada_data/hazards/Trial4_GB_*
%
%   copy single data to cluster:scp -r Documents/_GIT/climada_data/isimip/tc_tracks/Trial3_GB_dkgfdl_piControlcal dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/isimip/tc_tracks/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip04
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v04/.
%   copy results to dkrz:      scp -r /cluster/work/climate/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.
%
%   other option, a LSF pool, see http://www.clusterwiki.ethz.ch/brutus/Parallel_MATLAB_and_Brutus
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip04
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip04
%
%   run_on_desktop=1; % to test the job on a desktop
%   job_isimip04
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling job_isimip04 (all in
%       MATLAB command window), it sets the number of parallel pool workers
%       to two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.
%       Default=0 for cluster.
% OUTPUTS:
%   to disk, see PARAMETERS and climada folders
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20171025, copy from job_isimip03
% David N. Bresch, dbresch@ethz.ch, 20171027, scp added
% David N. Bresch, dbresch@ethz.ch, 20171029, new simulation names, job terminated 20171029_0908
% David N. Bresch, dbresch@ethz.ch, 20180312, redone for climada_global.tc.default_min_TimeStep=0.25; % 15 min
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
    %'Trial4_GB_dkgfdl_20thcal' % done 20180313, 12:53h
    'Trial4_GB_dkgfdl_piControlcal' % _gb_ to _GB_
    'Trial4_GB_dkipsl_20thcal'
    'Trial4_GB_dkipsl_piControlcal'
    'Trial4_GB_dkipsl_rcp26cal'
    'Trial4_GB_dkmiroc_20thcal'
    'Trial4_GB_dkmiroc_piControlcal'
    'Trial4_GB_dkmiroc_rcp26cal'
    % % in ascending size
    %     'Trial4_GB_dkgfdl_20thcal'
    %     'Trial4_GB_dkmiroc_20thcal'
    %     'Trial4_GB_dkipsl_20thcal'
    %     'Trial4_GB_dkmiroc_piControlcal'
    %     'Trial4_GB_dkgfdl_piControlcal' % _gb_ to _GB_
    %     'Trial4_GB_dkipsl_piControlcal'
    %     'Trial4_GB_dkmiroc_rcp26cal'
    %     'Trial4_GB_dkipsl_rcp26cal'
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
% fprintf('--> copy paste this into track_files in %s\n','job_isimip04');

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

% copy results to dkrz (no closing ; to log success):
% ---------------------
%[status,result]=system('scp -r    /cluster/scratch/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')
%[status,result] =system('scp -r -v /cluster/scratch/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.')

% %
% % SPECIAL CODE2 to inspect results
% % ----------------------------------
% result_hazard_dir='/Users/bresch/polybox/isimip/hazards_v04';
% fig_dir='/Users/bresch/Desktop/isimip';fig_ext='png';
% if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
% dd=dir([result_hazard_dir filesep '*.mat']);
% params.figure_scale=0; % no geographical scale on figure
% params.blue_ocean=1; % no geographical scale on figure
% for i=1:length(dd)
%     if ~dd(i).isdir && length(dd(i).name)>2
%         fprintf('processing %s: (%i of %i)\n',dd(i).name,i,length(dd))
%         hazard=climada_hazard_load([result_hazard_dir filesep dd(i).name]);
%         figure;climada_hazard_plot_nogrid(hazard,0,[],params); % plot max intensity
%         fN=strrep(dd(i).name,'.mat','');xlim([-180 180]);
%         ylim1=ylim;if ylim1(2)>90,ylim1(2)=90;end;if ylim1(2)<-90,ylim1(2)=-90;end
%         title_str=strrep(fN,'_',' ');ylim(ylim1);
%         title([title_str ' (max intensity)']);
%         saveas(gcf,[fig_dir filesep fN '.' fig_ext],fig_ext);
%         track_file=strrep(fN,        '_N_0360as','');
%         track_file=strrep(track_file,'_S_0360as','');
%         hemisphere='both';
%         if ~isempty(strfind(fN,'_S_')),hemisphere='S';end
%         if ~isempty(strfind(fN,'_N_')),hemisphere='N';end
%         tc_track=isimip_tc_track_load(track_file,hemisphere,180,-1);
%         figure;climada_tc_track_info(tc_track); % plot tracks
%         title(title_str);xlim([-180 180]);ylim(ylim1);
%         saveas(gcf,[fig_dir filesep track_file '_tc_track_' hemisphere '.' fig_ext],fig_ext);
%         close all
%     end % ~dd(i).isdir
% end % i

if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit