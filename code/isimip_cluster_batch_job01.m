% bath job for cluster, currently runs on Euler at ETH, see http://www.clusterwiki.ethz.ch/brutus/Getting_started_with_Euler 
% start with (if with parpool, 24 w): bsub -R "rusage[mem=5000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip01
cd /cluster/home/dbresch/climada
startup % climada_global exists afterwards
pwd % just to check where the job is running from
N_pool_workers=24; % for parpool
climada_global.parfor=1; % for parpool

% all tc_track files in ascending size
% 278954947 Mar 31 09:58 ../climada_data/isimip/temp_miroc20thcal.mat
% 296315458 Mar 31 09:52 ../climada_data/isimip/temp_ccsm420thcal.mat
% 307727612 Mar 31 09:53 ../climada_data/isimip/temp_gfdl520thcal.mat
% 311805989 Mar 31 09:54 ../climada_data/isimip/temp_hadgem20thcal.mat
% 316919214 Mar 31 10:01 ../climada_data/isimip/temp_mri20thcal.mat
% 323163248 Mar 31 09:59 ../climada_data/isimip/temp_mpi20thcal.mat
% 437561110 Mar 31 09:59 ../climada_data/isimip/temp_mirocrcp85cal_full.mat
% 458850472 Mar 31 09:52 ../climada_data/isimip/temp_ccsm4rcp85_full.mat
% 487672434 Mar 31 09:57 ../climada_data/isimip/temp_hadgemrcp85cal_full.mat
% 493729920 Mar 31 10:02 ../climada_data/isimip/temp_mrircp85cal_full.mat
% 494440909 Mar 31 09:54 ../climada_data/isimip/temp_gfdl5rcp85cal_full.mat
% 516923249 Mar 31 10:00 ../climada_data/isimip/temp_mpircp85cal_full.mat

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
%     'temp_ccsm420thcal' % all tc_track files alphabetical
%     'temp_mirocrcp85cal_full'
%     'temp_ccsm4rcp85_full'
%     'temp_mpi20thcal'
%     'temp_gfdl520thcal'
%     'temp_gfdl5rcp85cal_full'
%     'temp_mpircp85cal_full'
%     'temp_hadgem20thcal'
%     'temp_hadgemrcp85cal_full'
%     'temp_mri20thcal'
%     'temp_miroc20thcal'
%     'temp_mrircp85cal_full'
    };

% prepare centroids (full globe leads to segmentation fault)
centroids_S=climada_centroids_load('GLB_NatID_grid_0360as_adv_1');

centroids_N = centroids_S; % Northern Hemisphere
lat_pos=find(centroids_N.lat>0);
centroids_N.lon=centroids_N.lon(lat_pos);
centroids_N.lat=centroids_N.lat(lat_pos);
centroids_N.centroid_ID=centroids_N.centroid_ID(lat_pos);
centroids_N.distance2coast_km=centroids_N.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_N.NatID))=[];
centroids_N.NatID=centroids_N.NatID(lat_pos);

lat_pos=find(centroids_S.lat<=0); % Southern Hemisphere
centroids_S.lon=centroids_S.lon(lat_pos);
centroids_S.lat=centroids_S.lat(lat_pos);
centroids_S.centroid_ID=centroids_S.centroid_ID(lat_pos);
centroids_S.distance2coast_km=centroids_S.distance2coast_km(lat_pos);
lat_pos(lat_pos>length(centroids_S.NatID))=[];
centroids_S.NatID=centroids_S.NatID(lat_pos);

pool=parpool(N_pool_workers);
for file_i=1:length(track_files)
    
    tc_track=isimip_tc_track_load(track_files{file_i},'N',180,-1);
    hazard_name=[track_files{file_i} '_N_0360as'];
    isimip_tc_hazard_set(tc_track,hazard_name,centroids_N,0,hazard_name);
    
    tc_track=isimip_tc_track_load(track_files{file_i},'S',180,-1);
    hazard_name=[track_files{file_i} '_S_0360as'];
    isimip_tc_hazard_set(tc_track,hazard_name,centroids_S,0,hazard_name);
    
end % file_i
delete(pool)

exit
