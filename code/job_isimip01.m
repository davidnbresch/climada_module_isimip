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
    hazard_name=[filesep track_files{file_i} '_N_0360as'];
    hazard_set_file=[scratch_dir filesep hazard_name];
    isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_N,0,hazard_name);
    
    tc_track=isimip_tc_track_load(track_files{file_i},'S',180,-1); % Southern hemisphere
    hazard_name=[track_files{file_i} '_S_0360as'];
    hazard_set_file=[scratch_dir filesep hazard_name];
    isimip_tc_hazard_set(tc_track,hazard_set_file,centroids_S,0,hazard_name);
    
end % file_i
delete(pool)

exit % the cluster appreciates this, gives back memory etc.