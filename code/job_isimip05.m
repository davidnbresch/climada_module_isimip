% batch job for cluster: bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
% MODULE:
%   isimip
% NAME:
%   job_isimip05
% PURPOSE:
%   generate isimip tropical cyclone (TC) hazard event sets based on Kerry
%   Emmanuel TC track files i.e. files such as Trial4_GB_dkgfdl_20thcal.mat
%
%   But use a coarse (1deg) grid instead of the (default 360as one)
%
%   At the end, copy the hazard event sets to dkrz (isimip)
%
%   Once all runs are done, consider using SPECIAL CODE2 at the end below
%   to create check plots
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/job_isimip05.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy all data to cluster:  scp -r Documents/_GIT/climada_data/isimip/tc_tracks dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_data/isimip/.
%   copy single data to cluster:scp -r Documents/_GIT/climada_data/isimip/tc_tracks/Trial3_GB_dkgfdl_piControlcal dbresch@euler.ethz.ch:/cluster/home/dbresch/climada_data/isimip/tc_tracks/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v04/.
%   copy results to dkrz:      scp -r /cluster/scratch/dbresch/climada_data/hazards/*.mat b380587@mistralpp.dkrz.de:/work/bb0820/scratch/b380587/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip05
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20171123, copy from job_isimip04
%-

% PARAMETERS
%
FAST_TEST=0; % defalt=0, if =1, set -R "rusage[mem=500]"
%
% one should only have to edit this section
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada
scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
%scratch_dir = '/Users/bresch/Documents/_GIT/climada_data/isimip/scratch'; % for local tests
%
% the list of TC track files to be processed (see SPECIAL CODE below)
track_files={
    'Trial1_GB_dkmiroc_20thcal'
    'Trial1_GB_dkmiroc_rcp60cal'
    };

startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool

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

pool=parpool(N_pool_workers);
for file_i=1:length(track_files)
    
    hazard_name=[track_files{file_i} '_1deg_prob'];
    tc_track=isimip_tc_track_load(track_files{file_i},'both',180,-1); % both hemispheres
    
    % create probabilistic tracks
    [~,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,0); % wind speed decay at track nodes after landfall
    tc_track  = climada_tc_random_walk(tc_track); % overwrites tc_track to save memory
    tc_track  = climada_tc_track_wind_decay(tc_track, p_rel,0); % add the inland decay correction to all probabilistic nodes
    
    tc_track_prob_name=[scratch_dir filesep hazard_name '_tracks.mat'];
    fprintf('storing prob tracks %s\n',tc_track_prob_name)
    save(tc_track_prob_name,'tc_track',climada_global.save_file_version);
    
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
                        
    if FAST_TEST,tc_track=tc_track(1:1000);end % small subset for TEST
    fprintf('using %i centroids, %i (prob) tracks\n',length(centroids.lon),length(tc_track))

    hazard_set_file=[scratch_dir filesep hazard_name];
    isimip_tc_hazard_set(tc_track,hazard_set_file,centroids,0,hazard_name);
    
end % file_i
delete(pool)

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

exit % the cluster appreciates this, gives back memory etc.