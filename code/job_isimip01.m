% batch job for cluster: bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip01
% MODULE:
%   isimip
% NAME:
%   job_isimip01
% PURPOSE:
%   generate isimip tropical cyclone (TC) hazard event sets based on Kerry
%   Emmanuel TC track files
%
%   See PARAMETERS before running this
%
%   Once all runs are done, consider using SPECIAL CODE2 at the end below
%   to create check plots
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/euler_jobs/job_isimip01.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   run on cluster:            bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip01
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat Documents/_GIT/climada_data/hazards/.
%   copy results back polybox: scp -r dbresch@euler.ethz.ch:/cluster/scratch/dbresch/climada_data/hazards/*.mat /Users/bresch/polybox/isimip/hazards_v01/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip01
% EXAMPLE:
%   bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip01
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to scratch disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20170403
% David N. Bresch, dbresch@ethz.ch, 20170704, just checked
% David N. Bresch, dbresch@ethz.ch, 20170705, cluster hints added
% David N. Bresch, dbresch@ethz.ch, 20170827, last executed
%-

% PARAMETERS
% one should only have to edit this section
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada
scratch_dir = '/cluster/scratch/dbresch/climada_data/hazards';
%
% the list of TC track files to be processed
track_files={
    'temp_miroc20thcal' % all tc_track files in ascending size
    'temp_ccsm420thcal'
    'temp_gfdl520thcal'
    'temp_hadgem20thcal'
    'temp_mri20thcal'
    'temp_mpi20thcal'
    'temp_mirocrcp85cal_full'
    'temp_ccsm4rcp85_full'
    'temp_hadgemrcp85cal_full'
    'temp_mrircp85cal_full'
    'temp_gfdl5rcp85cal_full'
    'temp_mpircp85cal_full'
    };

startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool

% perpare centroids (full globe leads to segmentation fault)
centroids_S=climada_centroids_load('GLB_NatID_grid_0360as_adv_1');

centroids_N = centroids_S; % Northern hemisphere
lat_pos=find(centroids_N.lat>0);
centroids_N.lon=centroids_N.lon(lat_pos);
centroids_N.lat=centroids_N.lat(lat_pos);
centroids_N.centroid_ID=centroids_N.centroid_ID(lat_pos);
centroids_N.distance2coast_km=centroids_N.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_N.NatID))=[];
centroids_N.NatID=centroids_N.NatID(lat_pos);

lat_pos=find(centroids_S.lat<=0); % Southern hemisphere
centroids_S.lon=centroids_S.lon(lat_pos);
centroids_S.lat=centroids_S.lat(lat_pos);
centroids_S.centroid_ID=centroids_S.centroid_ID(lat_pos);
centroids_S.distance2coast_km=centroids_S.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_S.NatID))=[];
centroids_S.NatID=centroids_S.NatID(lat_pos);

pool=parpool(N_pool_workers);
for file_i=1:length(track_files)
    
    tc_track=isimip_tc_track_load(track_files{file_i},'N',180,-1); % Northern hemisphere
    hazard_name=[track_files{file_i} '_N_0360as'];
    hazard_set_file=[scratch_dir filesep hazard_name];
    isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_N,0,hazard_name);
    
    tc_track=isimip_tc_track_load(track_files{file_i},'S',180,-1); % Southern hemisphere
    hazard_name=[track_files{file_i} '_S_0360as'];
    hazard_set_file=[scratch_dir filesep hazard_name];
    isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_S,0,hazard_name);
    
end % file_i
delete(pool)

%
% % SPECIAL CODE2 to inspect results
% % ----------------------------------
% result_hazard_dir='/Users/bresch/polybox/isimip/hazards_v01';
% fig_dir='/Users/bresch/Desktop/isimip';fig_ext='png';
% if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
% dd=dir([result_hazard_dir filesep 'temp*.mat']);
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