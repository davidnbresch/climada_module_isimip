% isimip climada tests
% MODULE:
%   isimip
% NAME:
%   isimip_windfield_TESTS
% PURPOSE:
%   a batch job which runs several TESTS to compare windfields (and parts thereof)
%
%   One might rather run specifc segments (see code). Each segment
%   (starting after a line of =====) is independent of the code above it
%
%   see isimip_windfield_holland and climada_tc_windfield
%   next: isimip_tc_hazard_set
% CALLING SEQUENCE:
%   isimip_windfield_TESTS
% EXAMPLE:
%   isimip_windfield_TESTS
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   plots and stdout
% RESTRICTIONS:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161225, initial
% David N. Bresch, david.bresch@gmail.com, 20170202, climada_global.parfor
%-

global climada_global
climada_global.parfor=1;

% PARAMETERS
%
% easy way to switch on/off segments of tests - see below for what the
% switches mean
check_vtrans=0;
check_haversine=0;
check_Andrew=1;
check_historic_hazard_set=0;
check_probabilistic_hazard_set=0;

% check the translation windspeed (vtrans) calculation
% ====================================================
if check_vtrans
    tc_track=climada_subarray(isimip_ibtracs_read('TEST'),30:33);
    TEST_lon=-83-1:0.25:-79+1;TEST_lat= [24-1:0.25:28+1]; % Florida
    [X,Y] = meshgrid(TEST_lon,TEST_lat); % get 2-D
    centroids.lon=reshape(X,[1 numel(X)]); % make 1-D
    centroids.lat=reshape(Y,[1 numel(Y)]); % make 1-D
    
    %construct a circular track to check all
    p=-pi:pi/12:pi;
    tc_track.lon=2*sin(p)+mean(TEST_lon);
    tc_track.lat=2*cos(p)+mean(TEST_lat);
    fprintf('Note: you can run this for i=1..%i\n',floor(length(tc_track.lon)/2-1));
    
    figure('Name','isimip_windfield_TESTS: vtrans','Position',[10 854 1587 482])
    
    for ii=1:2
        if ii==1 % simple way  to show to legs of the track
            i=3;
        else
            i=mod(i+15,length(tc_track.lon)-1);
        end
        
        node_dx=tc_track.lon(i+1)-tc_track.lon(i); % forward vector
        node_dy=tc_track.lat(i+1)-tc_track.lat(i);
        node_len=sqrt(node_dx^2+node_dy^2); % length of forward vector
        
        % rotate 90 degrees, i.e.
        % x2=x* cos(a)+y*sin(a), with a=pi/2,cos(a)=0,sin(a)=1
        % y2=x*-sin(a)+Y*cos(a), therefore
        node_tmp=node_dx;node_dx=node_dy;node_dy=-node_tmp;
        
        centroids_dlon=centroids.lon-tc_track.lon(i); % vector from center
        centroids_dlat=centroids.lat-tc_track.lat(i);
        centroids_len=sqrt(centroids_dlon.^2+centroids_dlat.^2); % length
        
        % scalar product, a*b=|a|*|b|*cos(phi), phi angle between vectors
        cos_phi=(centroids_dlon*node_dx+centroids_dlat*node_dy)./centroids_len/node_len;
        if tc_track.lat(i)<0;cos_phi=-cos_phi;end % southern hemisphere
        
        cos_phi_grid = griddata(centroids.lon,centroids.lat,cos_phi,X,Y); % interpolate to grid
        
        subplot(1,2,ii)
        pcolor(X,Y,cos_phi_grid);hold on;shading flat;axis equal;
        caxis([-1 1]);colormap(jet);colorbar;
        climada_plot_world_borders(1);title('cos(phi)');
        plot(tc_track.lon,tc_track.lat,'-k');
        plot(tc_track.lon,tc_track.lat,'xk');
        plot([tc_track.lon(i),tc_track.lon(i+1)],[tc_track.lat(i),tc_track.lat(i+1)],'-r');
        plot(tc_track.lon(i),tc_track.lat(i),'xr');
        plot([tc_track.lon(i),tc_track.lon(i)+node_dx],[tc_track.lat(i),tc_track.lat(i)+node_dy],':r');
        xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
        
    end % ii
end % check_vtrans

% check haversine versus faster distance calculation
% ==================================================
if check_haversine
    TEST_lon=-83-1:0.25:-79+1;TEST_lat=24-1:0.25:28+1; % Florida
    [X,Y] = meshgrid(TEST_lon,TEST_lat); % get 2-D
    centroids.lon=reshape(X,[1 numel(X)]); % make 1-D
    centroids.lat=reshape(Y,[1 numel(Y)]); % make 1-D
    cos_centroids_lat = cos(centroids.lat/180*pi); % calculate once for speedup
    
    % calculate distance to all centroids
    r_arr_haversine=isimip_haversine(centroids.lon,centroids.lat,mean(TEST_lon),mean(TEST_lat),0.001); % km
    % faster version, 8km off for distances of >1000km
    r_arr_fast = sqrt(...
        ((centroids.lon-mean(TEST_lon)).*cos_centroids_lat).^2+...
        ( centroids.lat-mean(TEST_lat)                    ).^2)*111.1; % approx. conversion into km
    
    figure('Name','isimip_windfield_TESTS: distance comparison','Position',[10 854 1587 482])
    subplot(1,3,1)
    r_arr_haversine_grid = griddata(centroids.lon,centroids.lat,r_arr_haversine,X,Y); % interpolate back to grid
    pcolor(X,Y,r_arr_haversine_grid);hold on;shading flat;axis equal;
    climada_plot_world_borders(1);title('haversine (km)');
    plot(mean(TEST_lon),mean(TEST_lat),'ok','Markersize',10);
    plot(mean(TEST_lon),mean(TEST_lat),'xk','Markersize',10);
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
    
    subplot(1,3,2)
    r_arr_fast_grid = griddata(centroids.lon,centroids.lat,r_arr_fast,X,Y); % interpolate back to grid
    pcolor(X,Y,r_arr_fast_grid);hold on;shading flat;axis equal;
    climada_plot_world_borders(1);title('fast (km)');
    plot(mean(TEST_lon),mean(TEST_lat),'ok','Markersize',10);
    plot(mean(TEST_lon),mean(TEST_lat),'xk','Markersize',10);
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
    
    subplot(1,3,3)
    r_arr_fast_grid = griddata(centroids.lon,centroids.lat,r_arr_fast,X,Y); % interpolate back to grid
    pcolor(X,Y,r_arr_haversine_grid-r_arr_fast_grid);hold on;shading flat;axis equal;
    climada_plot_world_borders(1);title('haversine-fast (km)');
    plot(mean(TEST_lon),mean(TEST_lat),'ok','Markersize',10);
    plot(mean(TEST_lon),mean(TEST_lat),'xk','Markersize',10);
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
end % check_haversine

% compare windfield of Andrew (operatinal mode of wind field routines)
% ====================================================================
if check_Andrew
    % get tracks
    tc_track=climada_tc_track_load('tracks.atl_hist');
    tc_track=climada_tc_equal_timestep(tc_track,1/3); % 20 min
    % define centroids
    TEST_lon=-83-1:0.25:-79+1;TEST_lat=24-1:0.25:28+1; % Florida
    [X,Y] = meshgrid(TEST_lon,TEST_lat); % get 2-D
    centroids.lon=reshape(X,[1 numel(X)]);centroids.lat=reshape(Y,[1 numel(Y)]); % make 1-D
    t0=clock;gust_climada=climada_tc_windfield(    tc_track(1170),centroids);
    climada_time_sec=etime(clock,t0);
    fprintf('climada_tc_windfield     %2.4f sec\n',climada_time_sec);
    t0=clock;gust_isimip =isimip_windfield_holland(tc_track(1170),centroids);
    isimip_time_sec=etime(clock,t0);
    fprintf('isimip_windfield_holland %2.4f sec (factor %2.2f compared to climada)\n',isimip_time_sec,isimip_time_sec/climada_time_sec);
    
    figure('Name','isimip_windfield_TESTS: Andrew comparison','Position',[10 854 1587 482])
    
    subplot(1,3,1);
    gust_climada_grid = griddata(centroids.lon,centroids.lat,gust_climada,X,Y); % interpolate back to grid
    gust_climada_grid = gust_climada_grid/1.27;
    fprintf('NOTE: climada wind divided by gust factor 1.27 to obtain gust, not peak gust (as used in climada historically)\n');
    pcolor(X,Y,gust_climada_grid);hold on;shading flat;axis equal;caxis([0 100]);
    climada_plot_world_borders(1);title('climada (m/s)');
    plot(tc_track(1170).lon,tc_track(1170).lat,'.k')
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
    
    subplot(1,3,2);
    gust_isimip_grid = griddata(centroids.lon,centroids.lat,gust_isimip,X,Y); % interpolate back to grid
    pcolor(X,Y,gust_isimip_grid);hold on;shading flat;axis equal;caxis([0 100]);
    climada_plot_world_borders(1);title('isimip (m/s)');
    plot(tc_track(1170).lon,tc_track(1170).lat,'.k')
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
    
    subplot(1,3,3);
    pcolor(X,Y,gust_climada_grid-gust_isimip_grid);hold on;shading flat;axis equal;caxis([-20 20]);
    climada_plot_world_borders(1);title('climada-isimip (m/s)');
    plot(tc_track(1170).lon,tc_track(1170).lat,'.k')
    xlim([min(TEST_lon),max(TEST_lon)]);ylim([min(TEST_lat),max(TEST_lat)]);colormap(jet);colorbar;
    
end % check_Andrew

% calculate historic haazrd event set
% ===================================
if check_historic_hazard_set
    
    % get tracks
    tc_track=climada_tc_track_load('tracks.atl_hist');
    tc_track=climada_tc_equal_timestep(tc_track,1/3); % 20 min
    
    entity=climada_entity_load('USA_UnitedStates_Florida');
    additional_lat=24:.5:32;
    additional_lon=[-79 -79.5];
    [X,Y] = meshgrid(additional_lon,additional_lat);
    additional_lon=reshape(X,[1 numel(X)]);additional_lat=reshape(Y,[1 numel(Y)]); % make 1-D
    centroids.lon=[entity.assets.lon additional_lon];
    centroids.lat=[entity.assets.lat additional_lat];
    centroids.centroid_ID=1:numel(centroids.lon);
    plot(centroids.lon,centroids.lat,'.r');hold on;
    climada_plot_world_borders(1);title('centroids');
    xlim([min(centroids.lon),max(centroids.lon)]);ylim([min(centroids.lat),max(centroids.lat)]);
    
    hazard_isimip_hist  = isimip_tc_hazard_set( tc_track,'_TC_atl_isimip_hist', centroids);
    hazard_climada_hist = climada_tc_hazard_set(tc_track,'_TC_atl_climada_hist',centroids);
    hazard_climada_hist.intensity=hazard_climada_hist.intensity/1.27; % gust, not peak gust to compare with isimip
    
    hazard_climada_hist.intensity(hazard_climada_hist.intensity<20)=0;
    hazard_isimip_hist.intensity( hazard_isimip_hist.intensity <20)=0;
    
    res_isimip_hist =climada_hazard_check(hazard_isimip_hist);
    res_climada_hist=climada_hazard_check(hazard_climada_hist,res_isimip_hist);
    
    figure('Name','isimip_windfield_TESTS: hazard stats (climada, historic)','Position',[20 268 1481 768])
    climada_hazard_stats(hazard_climada_hist);
    
    figure('Name','isimip_windfield_TESTS: hazard stats (isimip, historic)','Position',[20 268 1481 768])
    climada_hazard_stats(hazard_isimip_hist);
    
end % check_historic_hazard_set

% calculate probabilistic haazrd event set
% ========================================
if check_probabilistic_hazard_set
    
    % get tracks
    tc_track=climada_tc_track_load('tracks.atl_hist');
    
    % construct probabilistic ones
    [~,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,0); % pressure decay relation
    tc_track  = climada_tc_random_walk(tc_track); % overwrites tc_track to save memory
    tc_track  = climada_tc_track_wind_decay(tc_track, p_rel,0);
    
    tc_track=climada_tc_equal_timestep(tc_track,1/3); % 20 min
    
    entity=climada_entity_load('USA_UnitedStates_Florida');
    additional_lat=24:.5:32;
    additional_lon=[-79 -79.5];
    [X,Y] = meshgrid(additional_lon,additional_lat);
    additional_lon=reshape(X,[1 numel(X)]);additional_lat=reshape(Y,[1 numel(Y)]); % make 1-D
    centroids.lon=[entity.assets.lon additional_lon];
    centroids.lat=[entity.assets.lat additional_lat];
    centroids.centroid_ID=1:numel(centroids.lon);
    plot(centroids.lon,centroids.lat,'.r');hold on;
    climada_plot_world_borders(1);title('centroids');
    xlim([min(centroids.lon),max(centroids.lon)]);ylim([min(centroids.lat),max(centroids.lat)]);
    
    hazard_isimip_prob  = isimip_tc_hazard_set( tc_track,'_TC_atl_isimip_prob', centroids);
    hazard_climada_prob = climada_tc_hazard_set(tc_track,'_TC_atl_climada_prob',centroids);
    
    % NOTE: climada scaled, NOT in saved hazard set
    hazard_climada_prob.intensity=hazard_climada_prob.intensity/1.27; % gust, not peak gust to compare with isimip
    
    hazard_climada_prob.intensity(hazard_climada_prob.intensity<20)=0;
    hazard_isimip_prob.intensity(hazard_isimip_prob.intensity<20)=0;
    
    res_isimip_prob =climada_hazard_check(hazard_isimip_prob);
    res_climada_prob=climada_hazard_check(hazard_climada_prob,res_isimip_prob);
    
    figure('Name','isimip_windfield_TESTS: hazard stats (climada, probabilistic)','Position',[20 268 1481 768])
    climada_hazard_stats(hazard_climada_prob);set(gcf,'Position',[20 268 1481 768]);
    
    figure('Name','isimip_windfield_TESTS: hazard stats (isimip, probabilistic)','Position',[20 268 1481 768])
    climada_hazard_stats(hazard_isimip_prob);set(gcf,'Position',[20 268 1481 768]);
    
end % check_probabilistic_hazard_set



tc_track_hist=climada_tc_read_unisys_database('nio'); % nio: North Indian Ocean
climada_tc_track_info(tc_track_hist,1);xlim([40 110]);ylim([0 40]) % Fig 1
tc_track_prob9=climada_tc_random_walk(tc_track_hist,9);
climada_tc_track_info(tc_track_prob9,1); xlim([40 110]);ylim([0 40]) % Fig 2
p.resolution_km=10;entity=climada_nightlight_entity('India','Maharashtra',p);
entity.assets.Value=entity.assets.Value/sum(entity.assets.Value)*250e9*(2+1);entity.assets.Cover=entity.assets.Value;
climada_entity_plot(entity,3) % Fig 3
hazard_hist=climada_tc_hazard_set(tc_track_hist,'NOSAVE',entity,1);
climada_hazard_stats(hazard_hist) % Fig 4
hazard_prob9=climada_tc_hazard_set(tc_track_prob9,'NOSAVE',entity,1);
climada_hazard_stats(hazard_prob9) % Fig 5
tc_track_prob99=climada_tc_random_walk(tc_track_hist,99);
hazard_prob99=climada_tc_hazard_set(tc_track_prob99,'NOSAVE',entity,1);
climada_hazard_stats(hazard_prob99) % Fig 6

