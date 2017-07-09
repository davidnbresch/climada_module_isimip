function [gust,hazard]=isimip_windfield_holland(tc_track,centroids,~,silent_mode,~)
% TC windfield calculation
% MODULE:
%   isimip
% NAME:
%   isimip_windfield_holland
% PURPOSE:
%   given a TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the wind field at locations (=centroids). 
%
%   Using a FAST distance calculation (not haversine), which speeds
%   calculations up by factor 2-3. See r_arr and dist in code.
%
%   see isimip_windfield_TESTS to check for haversine/fast difference, for
%   the vtrans calculation and for the difference between isimip and
%   default climada windfield calculation
%
%   see TEST for testing options (to check for 1-to-1 correspondence with
%   python). There are several TEST sections, indicated always with
%   TEST_START and ending with TEST_END
%
%   units not double-checked yet, see UNIT_TEST
%
%   See also the original (not fully functional) Matlab version of
%   helper_advanced_windfield_global.py, i.e.
%   helper_advanced_windfield_global.m
%
%   so far adapted to windfield by holland1980, holland2008, holland2010 cf peduzzi2012
%   can be easily expanded to different wind-fields eg chavez2015.
%   Greg J. Holland, 1980: An Analytic Model of the Wind and Pressure Profiles in
%   Hurricanes. DOI: http://dx.doi.org/10.1175/1520-0493(1980)108<1212:AAMOTW>2.0.CO;2 
%
% CALLING SEQUENCE:
%   gust = isimip_windfield_holland(tc_track,centroids,~,silent_mode,~)%
% EXAMPLE:
%   gust=isimip_windfield_holland; % TEST, does the same as next line
%   gust=isimip_windfield_holland(climada_subarray(isimip_ibtracs_read('TEST'),30:33));
%   tc_track=isimip_ibtracs_read('TEST');
%   tc_track=climada_tc_equal_timestep(tc_track,1); % make equal timestep, 1 hour
%   gust=isimip_windfield_holland(tc_track)
%   climada_color_plot(gust,centroids.lon,centroids.lat);
%
%   a bit special to plot single event footprint
%   tc_track=climada_tc_track_load([climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat']);
%   for ii=1:length(tc_track),iii=findstr('2005236N23285',tc_track(ii).ID_str);if iii==1,i=ii;end;end % Katrina
%   %for ii=1:length(tc_track),iii=findstr('2013306N07162',tc_track(ii).ID_str);if iii==1,i=ii;end;end % Haiyan
%   [gust,hazard] = isimip_windfield_holland(climada_tc_equal_timestep(tc_track(i)),'0360as');
%   climada_hazard_plot_nogrid(hazard);title('2005236N23285');
% INPUTS:
%   tc_track: a structure with the single track information (length(tc_track)!=1)
%       see e.g. climada_tc_read_unisys_tc_track
%       tc_track.Azimuth and/or tc_track.Celerity calculated, if not existing
%       but climada_tc_equal_timestep mist have been run and
%       tc_track.MaxSustainedWind must exist on input
%       If empty: TEST mode, use three nodes of Andrew
%   centroids: a structure with the centroids information (see e.g.
%       climada_centroids_read):
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
%       If empty together with tc_track: TEST mode, use three nodes of Andrew
%       If ='0360as' or ='0150as', generate centroids to cover the track
%       region (mainly used for nice footprint plots, also over water)
% OPTIONAL INPUT PARAMETERS:
%   silent_mode: default=1 (for speedup), do not print or plot
%       if =-7, do show output to stdout and print windfield, only makes
%       sense for TEST mode, since comparison windfields are read from files
% OUTPUTS:
%   gust: the 1-im peak gust windfield [m/s] at all centroids, NOT sparse 
%       for speedup, i.e. convert like hazard.intensity()=sparse(gust)
%   hazard: the single-event hazard set, if requested (not provied on
%       output by default for speedup)
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Tobias Geiger, geiger@pik-potsdam.de, 2016, Copyright, python verison
% David N. Bresch, david.bresch@gmail.com, 20161204, intial
% David N. Bresch, david.bresch@gmail.com, 20161205, traslation added
% David N. Bresch, david.bresch@gmail.com, 20161223, traslation removed for TESTS
% David N. Bresch, david.bresch@gmail.com, 20161224, TEST options finalized
% David N. Bresch, david.bresch@gmail.com, 20161225, final output in m/s
% David N. Bresch, david.bresch@gmail.com, 20170707, centroids='0360as' and hazard as optional  output
%-

if ~exist('tc_track' ,'var')
    tc_track=climada_subarray(isimip_ibtracs_read('TEST'),30:33);
    fprintf('WARNING: simple TEST MODE in %s (using a few nodes of Andrew)\n',mfilename);
end
if ~exist('centroids','var'),centroids=[]; end
if ~exist('silent_mode','var'),silent_mode=1; end

% PARAMETERS
%
% only deal with centtroids not further away than from any node
centroid_node_max_dist_km=300;
%
% minimum windspeed stored (for storage management reasons, i.e. as later
% used as sparse matrix, this saves huge amounts of memory)
min_wind_threshold=34; % knots (UNIT_TEST)
%
% the model to be used
model='H08';
%
rho=1.15; % data for windfield calculation
EnvironmentalPressure_default=1010; % mb
%
% extrapolate RadiusMaxWind from pressure if not given
rmax1=15;rmax2=25;rmax3=50; % rmax-thresholds in nm
p1=950;p2=980;p3=1020; % in mb
%
% next line the test area (Southern tip of Florida)
% no need ot comment out if no test
test_lon=-83:0.1:-79;test_lat=24:0.1:28;
%
% the degrees we add 'around' the track in case we generate centroids
% (case centroids='0360as'...)
d=3; % default=3
%
% % TEST_START
% global climada_global
% if ~climada_init_vars,return;end % init/import global variables
% % here a section reading TEST output from Python starts
% TEST_comparison_file=[climada_global.data_dir filesep 'isimip' ...
%     filesep 'ibtracs_andrew_test' filesep 'ibtracs_andrew_test_3nodes_with-vtrans.nc'];
% %    filesep 'ibtracs_andrew_test' filesep 'ibtracs_andrew_test_3nodes_no-vtrans.nc'];
% TEST_comparison_file_nvt=[climada_global.data_dir filesep 'isimip' ...
%     filesep 'ibtracs_andrew_test' filesep 'ibtracs_andrew_test_3nodes_no-vtrans.nc'];
% if exist('TEST_comparison_file','var')
%     fprintf('TEST mode, reading %s\n',TEST_comparison_file);
%     nc.vmax=ncread(TEST_comparison_file,'vmax');
%     nc.vmax_nvt=ncread(TEST_comparison_file_nvt,'vmax');
%     nc.vtrans=nc.vmax-nc.vmax_nvt;
%     nc.lon=ncread(TEST_comparison_file,'lon');
%     nc.lat=ncread(TEST_comparison_file,'lat');
%     [nc.gridlon,nc.gridlat]=meshgrid(nc.lon,nc.lat);
%     nc.gridlon=nc.gridlon';nc.gridlat=nc.gridlat';
%     test_lon=nc.lon;
%     test_lat=nc.lat;
%     %XLim=[-83 -79];YLim=[24 28];
%     %XLim=[-88 -75];YLim=[20 32];
%     XLim=[-82 -77];YLim=[23 28]; % the area of the two nodes we inspect
%     silent_mode=-7; % force output to stdout and plots
% end
% % TEST_END

if isempty(centroids) % generate test centroids
    [centroids.lon,centroids.lat] = meshgrid(test_lon,test_lat); % get 2-D
    centroids.lon=reshape(centroids.lon,[1 numel(centroids.lon)]); % make 1-D
    centroids.lat=reshape(centroids.lat,[1 numel(centroids.lat)]); % make 1-D
end

if ischar(centroids)
    if strcmp(centroids,'0360as')
        ddd=0.025;
    elseif strcmp(centroids,'0150as')
        ddd=0.0125;
    else
        fprintf('centroids=%s not implemented\n',centroids);
        return
    end
    centroids=[];
    xgv=min(tc_track.lon)-d:ddd:max(tc_track.lon)+d;
    ygv=min(tc_track.lat)-d:ddd:max(tc_track.lat)+d;
    [X,Y] = meshgrid(xgv,ygv);
    centroids.lon=reshape(X,[1 numel(X)]); % make 1-D
    centroids.lat=reshape(Y,[1 numel(Y)]); % make 1-D
    
end

% use default EnvironmentalPressure, if not provided
if ~isfield(tc_track,'EnvironmentalPressure')
    tc_track.EnvironmentalPressure=tc_track.CentralPressure*0+EnvironmentalPressure_default;
end

if ~isfield(tc_track,'RadiusMaxWind') % set with default
    tc_track.RadiusMaxWind=tc_track.CentralPressure*0;
end

t0=clock;

% make sure that CentralPressure never exceeds EnvironmentalPressure
% (original: if pcen > penv,pcen=penv;end with pcen=tc_track.CentralPressure)
p_points=tc_track.CentralPressure>tc_track.EnvironmentalPressure;
tc_track.CentralPressure(p_points)=tc_track.EnvironmentalPressure(p_points);

% extrapolate RadiusMaxWind from pressure if not given
% (original: see extra_rmax)
% cubic fit , ibtracs 1980-2013 (r2=0.22)
% old: rmax_from_p= -45201.105207 + 146.043463*pcen - 0.157263*pcen**2 + 0.000056*pcen**3 + 0.097515*lat + 0.016056*lon
% rmax extrapolated from historic distribution, rmax-thresholds in nm
tc_track.RadiusMaxWind(tc_track.CentralPressure<=p1)=rmax1;
p_points=tc_track.CentralPressure>p1 & tc_track.CentralPressure<=p2;
tc_track.RadiusMaxWind(p_points)=(tc_track.CentralPressure(p_points)-p1)*(rmax2-rmax1)/(p2-p1) + rmax1;
p_points=tc_track.CentralPressure>p2;
tc_track.RadiusMaxWind(p_points)=(tc_track.CentralPressure(p_points)-p2)*(rmax3-rmax2)/(p3-p2) + rmax2;
tc_track.RadiusMaxWind=tc_track.RadiusMaxWind*1.852; % nautical mile to km

n_nodes=length(tc_track.lon);

% add one more point to the end of the track (for forward vector)
tc_track.lon(end+1)=tc_track.lon(end);
tc_track.lat(end+1)=tc_track.lat(end);

gust=centroids.lon*0; % init output

cos_centroids_lat = cos(centroids.lat/180*pi); % calculate once for speedup

for node_i=2:n_nodes % process each node (later speedup potential)
    
    % calculate distance to all centroids
    %r_arr=isimip_haversine(centroids.lon,centroids.lat,tc_track.lon(node_i),tc_track.lat(node_i),0.001); % km
    % faster version, 8km off for distances of >1000km
    r_arr = sqrt(...
        ((centroids.lon-tc_track.lon(node_i)).*cos_centroids_lat).^2+...
        ( centroids.lat-tc_track.lat(node_i)                    ).^2)*111.1; % approx. conversion into km
       
    %cc=true(1,length(centroids.lon)); % use all, cc stands for Centroids Close enough
    cc=find(r_arr<centroid_node_max_dist_km); % only centroids close enough
    r_arr=r_arr(cc); % restrict r_arr
    
    % calculate translation speed of track
    %dist=isimip_haversine(tc_track.lon(node_i-1),tc_track.lat(node_i-1),tc_track.lon(node_i),tc_track.lat(node_i));
    % faster version, difference to haversine does not matter for such distances
    dist=sqrt(((tc_track.lon(node_i-1)-tc_track.lon(node_i))*cos(tc_track.lat(node_i)/180*pi))^2+...
        (tc_track.lat(node_i-1)-tc_track.lat(node_i))^2)*111.1; % approx. conversion into km
    dist= dist/1.852; % dist to nautical miles
    % hours = hours between track coordiantes
    vtrans=dist/tc_track.TimeStep(node_i); % nautical miles/hour
    if vtrans>30,vtrans=30;end % limit to 30 nmph
    
    % convert variable names to the ones used in isimip (python)
    vmax    =tc_track.MaxSustainedWind(node_i);
    rmax    =tc_track.RadiusMaxWind(node_i);
    pcen    =tc_track.CentralPressure(node_i);
    prepcen =tc_track.CentralPressure(node_i-1); % previous node
    penv    =tc_track.EnvironmentalPressure(node_i);
    tint    =tc_track.TimeStep(node_i);
    %xc      =tc_track.lon(node_i); % not used
    yc      =tc_track.lat(node_i);
    
    % adjust pressure at previous track point
    if prepcen < 850,prepcen=pcen;end
    
    % calculate the symmetric Holland windfield
    % -----------------------------------------
    
    b=b_value(vtrans,vmax,penv,pcen,rho);
    xx=x_value(penv,pcen);
    bs=bs_value(vtrans,penv,pcen,prepcen,yc,xx,tint);
    
    if strcmp(model,'H80')
        vwind=stat_holland(r_arr,rmax,b ,penv,pcen,yc,rho);
    elseif strcmp(model,'H08')
        % same as H80 but with bs instead of b, additional info required previous track point central pressure for pressure gradient
        vwind=stat_holland(r_arr,rmax,bs,penv,pcen,yc,rho);
        %     elseif strcmp(model,'H08_v')
        %         % same as H08 but uses max-wind speed instead of pressure
        %         vwind=stat_holland_alt(r_arr,rmax,vmax,bs,yc);
        %     elseif strcmp(model,'H10')
        %         % in holland2010 x_exp varies with distance, not implemented
        %         % so far identical to holland2008 but without gradient wind
        %         % note: NO gradient wind implemented for holland2010!!!
        %         x_exp=0.5;
        %         vwind=stat_holland2010(r_arr,rmax,penv,pcen,bs,x_exp,rho);
    end
        
    % calculate angular field to add translational wind
    % -------------------------------------------------
    
    % figure which side of track, hence add/subtract translational wind
    node_dx=tc_track.lon(node_i+1)-tc_track.lon(node_i); % track forward vector
    node_dy=tc_track.lat(node_i+1)-tc_track.lat(node_i);
    node_len=sqrt(node_dx^2+node_dy^2); % length of track forward vector
    
    % we use the scalar product of the track forward vector and the vector
    % towards each centroid to figure the angle between and hence whether
    % the translational wind needs to be added (on the right side of the
    % track for Northern hemisphere) and to which extent (100% exactly 90
    % to the right of the track, zero in front of the track)
    
    % hence, rotate track forward vector 90 degrees clockwise, i.e.
    % x2=x* cos(a)+y*sin(a), with a=pi/2,cos(a)=0,sin(a)=1
    % y2=x*-sin(a)+Y*cos(a), therefore
    node_tmp=node_dx;node_dx=node_dy;node_dy=-node_tmp;
    
    % the vector towards each centroid
    centroids_dlon=centroids.lon(cc)-tc_track.lon(node_i); % vector from center
    centroids_dlat=centroids.lat(cc)-tc_track.lat(node_i);
    centroids_len=sqrt(centroids_dlon.^2+centroids_dlat.^2); % length
    
    % scalar product, a*b=|a|*|b|*cos(phi), phi angle between vectors
    cos_phi=(centroids_dlon*node_dx+centroids_dlat*node_dy)./centroids_len/node_len;
    if tc_track.lat(node_i)<0;cos_phi=-cos_phi;end % southern hemisphere
    
    % calculate vtrans wind field array assuming that
    % - effect of vtrans decreases with distance from eye (r_arr_normed)
    % - vtrans is added 100% to the right of the track, 0% in front etc. (cos_phi)
    r_arr_normed=rmax./r_arr;
    r_arr_normed(r_arr_normed>1)=1;
    vtrans_arr=vtrans*r_arr_normed.*cos_phi;
    
    % add translational wind to obtain final windfield
    % ------------------------------------------------
    %vfull=vwind; % TEST: symmetric windfield only
    vfull=vwind+vtrans_arr; % symmetric windfield plus translation
    
    vfull(vfull<min_wind_threshold)=0;
    
    gust(cc)=max(gust(cc),vfull); % keep maximum instantaneous wind
    
end % node_i

if silent_mode==-7
    
    isimip_etime=etime(clock,t0);
    fprintf('isimip windfield took %f sec\n',isimip_etime);
    fprintf('plotting ...\n');
    
    % % SWITCH this section on to compare with and without vtrans
    % % run code FIRST with
    % % line about 246 vfull=vwind;             instead of
    % % line about 246 vfull=vwind+vtrans_arr;
    % % and uncomment ONLY the next line
    % %save('TEST_gust_nvt','gust'); return
    % % then siwthc back to vfull=vwind+vtrans_arr and uncomment code below
    % %
    % gust_nvt=load('TEST_gust_nvt'); % load gust with no vtrans
    % caxis_rng=[-10 10];XLim=[-82 -77];YLim=[23 28];
    % vtrans=gust-gust_nvt.gust;
    % [X,Y] = meshgrid(test_lon,test_lat); % get 2-D
    % gridvtrans = griddata(centroids.lon,centroids.lat,vtrans,X,Y); % interpolate to grid
    % figure('Name','vtrans comparison','Position',[10 854 1587 482])
    % subplot(1,3,1);
    % pcolor(X,Y,gridvtrans);hold on;shading flat;axis equal;
    % caxis(caxis_rng);climada_plot_world_borders(1);title('MATLAB');
    % plot(tc_track.lon,tc_track.lat,'-k');plot(tc_track.lon,tc_track.lat,'xk');
    % xlim(XLim);ylim(YLim);colormap(jet);colorbar;
    % subplot(1,3,2);
    % pcolor(nc.gridlon,nc.gridlat,nc.vtrans);hold on;shading flat;axis equal;
    % caxis(caxis_rng);climada_plot_world_borders(1);title('Python');
    % plot(tc_track.lon,tc_track.lat,'-k');plot(tc_track.lon,tc_track.lat,'xk');
    % xlim(XLim);ylim(YLim);colormap(jet);colorbar;
    % subplot(1,3,3) % the difference
    % d_vtrans=gridvtrans'-double(nc.vtrans);
    % pcolor(nc.gridlon,nc.gridlat,d_vtrans);hold on;shading flat;axis equal;
    % caxis(caxis_rng);climada_plot_world_borders(1);title('MATLAB-Python');
    % plot(tc_track.lon,tc_track.lat,'-k');plot(tc_track.lon,tc_track.lat,'xk');
    % xlim(XLim);ylim(YLim);colormap(jet);colorbar
    % return
    % % end of SWITCH this section on to compare with and without vtrans
    
    figure('Name','windfield comparison','Position',[10 854 1587 482])
    subplot(1,3,1)
    [X,Y] = meshgrid(test_lon,test_lat); % get 2-D
    gridgust = griddata(centroids.lon,centroids.lat,gust,X,Y); % interpolate to grid
    pcolor(X,Y,gridgust);hold on;shading flat;axis equal;
    caxis([0 140]);
    climada_plot_world_borders(1);title('MATLAB');
    plot(tc_track.lon,tc_track.lat,'-k');
    plot(tc_track.lon,tc_track.lat,'xk');
    xlim(XLim);
    ylim(YLim);
    colormap(jet);colorbar
    gridgust(gridgust<34)=0;
    pos=find(gridgust>0);
    xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(gridgust(pos))),max(max(gridgust(pos)))));
    
    if exist('nc','var')
        subplot(1,3,2)
        pcolor(nc.gridlon,nc.gridlat,nc.vmax);hold on;shading flat;axis equal;
        caxis([0 140]);
        climada_plot_world_borders(1);title('netCDF (Python)');
        plot(tc_track.lon,tc_track.lat,'-k');
        plot(tc_track.lon,tc_track.lat,'xk');
        xlim(XLim);
        ylim(YLim);
        colormap(jet);colorbar
        pos=find(nc.vmax>0);
        xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(nc.vmax(pos))),max(max(nc.vmax(pos)))));
        
        subplot(1,3,3) % the difference
        d_vmax=gridgust'-double(nc.vmax);
        pcolor(nc.gridlon,nc.gridlat,d_vmax);hold on;shading flat;axis equal;
        caxis([-20 20]);
        climada_plot_world_borders(1);title('MATLAB-Python');
        plot(tc_track.lon,tc_track.lat,'-k');
        plot(tc_track.lon,tc_track.lat,'xk');
        xlim(XLim);
        ylim(YLim);
        colormap(jet);colorbar
        pos=find(abs(d_vmax)>eps);
        xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(d_vmax(pos))),max(max(d_vmax(pos)))));
        
    end % exist('nc','var')
    
    t0=clock;
    climada_gust = climada_tc_windfield(tc_track,centroids,[],-1); % in m/s
    climada_etime=etime(clock,t0);
    fprintf('climada windfield took %f sec\n',climada_etime);
    fprintf('plotting ...\n');

    figure('Name','windfield comparison','Position',[10 854 1587 482])
    subplot(1,3,1)
    climada_gridgust = griddata(centroids.lon,centroids.lat,climada_gust,X,Y); % interpolate to grid
    pcolor(X,Y,climada_gridgust);hold on;shading flat;axis equal;
    caxis([0 100]);colormap(jet);colorbar
    climada_plot_world_borders(1);title('climada (m/s)');
    plot(tc_track.lon,tc_track.lat,'-k');
    plot(tc_track.lon,tc_track.lat,'xk');
    xlim(XLim);ylim(YLim);
    %climada_gridgust(climada_gridgust<34)=0;
    pos=find(climada_gridgust>0);
    xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(climada_gridgust(pos))),max(max(climada_gridgust(pos)))));
    
    gridgust=gridgust*1.852/3.6; % from nautical miles/hour to km/h to m/s UNIT_TEST
    
    subplot(1,3,2)
    pcolor(X,Y,gridgust);hold on;shading flat;axis equal;
    caxis([0 100]);
    climada_plot_world_borders(1);title('isimip (m/s)');
    plot(tc_track.lon,tc_track.lat,'-k');
    plot(tc_track.lon,tc_track.lat,'xk');
    xlim(XLim);
    ylim(YLim);
    colormap(jet);colorbar
    pos=find(gridgust>0);
    xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(gridgust(pos))),max(max(gridgust(pos)))));
    
    subplot(1,3,3) % the difference
    d_gust=climada_gridgust-gridgust;
    pcolor(X,Y,d_gust);hold on;shading flat;axis equal;
    caxis([-20 20]);
    climada_plot_world_borders(1);title('climada-isimip (m/s)');
    plot(tc_track.lon,tc_track.lat,'-k');
    plot(tc_track.lon,tc_track.lat,'xk');
    xlim(XLim);
    ylim(YLim);
    colormap(jet);colorbar
    pos=find(abs(d_gust)>eps);
    xlabel(sprintf('min/max: %2.1f/%2.1f\n',min(min(d_gust(pos))),max(max(d_gust(pos)))));
    
end % silent_mode==-7

gust=gust*1.852/3.6; % from nautical miles/hour to km/h to m/s UNIT_TEST

if nargout==2
    % compile a single-event hazard set
    hazard.peril_ID='TC';
    hazard.units='m/s';
    hazard.lon=centroids.lon;
    hazard.lat=centroids.lat;
    hazard.intensity=gust;
    hazard.fraction=spones(hazard.intensity); % fraction 100%
    hazard.centroid_ID=1:length(hazard.lon);
    hazard.date=datestr(now);
    hazard.filename=mfilename;
    hazard.reference_year=tc_track.yyyy(1);
    hazard.event_ID=1;
    hazard.event_count=1;
    hazard.orig_event_flag=1;
    hazard.orig_event_count=1;
    hazard.orig_years=tc_track.yyyy(1);
    hazard.frequency=1;
end

end % isimip_windfield_holland


% follow helper functions (the ones easy to convert from python ;-)
% -----------------------------------------------------------------

function b=b_value(vtrans,vmax,penv,pcen,rho)
if vmax>0
    b=rho*   exp(1)*(0.51444444444*(vmax-vtrans)) ^2/(100*(penv-pcen));
    % b=rho*np.exp(1)*(0.51444444444*(vmax-vtrans))**2/(100*(penv-pcen))
    if b < 1
        b=1.0;
    elseif b > 2.5
        b=2.5;
    else
        b=1.0;
    end
else
    b=1; % ADDED dnb
end
end % b_value

function xx=x_value(penv,pcen)
xx=0.6*(1.-(penv-pcen)/215);
end % x_value

function bs=bs_value(vtrans,penv,pcen,prepcen,lat,xx,tint)
vt_ms=vtrans*0.51444444444;
bs=-4.4e-5*(penv-pcen) ^2 + 0.01*(penv-pcen) + 0.03*(pcen-prepcen)/tint - 0.014*   abs(lat) + 0.15*vt_ms ^xx + 1.0;
% bs=-4.4e-5*(penv-pcen)**2 + 0.01*(penv-pcen) + 0.03*(pcen-prepcen)/tint - 0.014*np.abs(lat) + 0.15*vt_ms**xx + 1.0
end % bs_value

function v_arr=stat_holland(r_arr,rmax,b,penv,pcen,lat,rho)
% holland symmetric and static wind field according to Holland1980 and with b=bs according to Holland2008
lat=abs(lat);
f=2*0.0000729 * sin(lat*pi/180);
% units are m/s
v_arr=r_arr*0; % init
for i=1:length(r_arr)
    v_arr(i)=  sqrt(((100*b/rho*(rmax./r_arr(i))^b*(penv-pcen)*exp(-(rmax/r_arr(i))^b))^0.5)^2+(1000*0.5*r_arr(i)*f)^2)-0.5*1000*r_arr(i)*f;
    % v_arr=   sqrt(((100*b/rho*(rmax./r_arr).^b*(penv-pcen)*  exp(-(rmax./r_arr).^b)).^0.5).^2+(1000*0.5*r_arr*f) ^2) - 0.5*1000*r_arr*f;
    % v_arr=np.sqrt(((100*b/rho*(rmax/r_arr)**b*(penv-pcen)*np.exp(-(rmax/r_arr)**b))**0.5)**2+(1000*0.5*r_arr*f)**2) - 0.5*1000*r_arr*f
end % i
v_arr(isnan(v_arr))=0;
% translate to knots
v_arr=v_arr/0.51444444444;
end % stat_holland