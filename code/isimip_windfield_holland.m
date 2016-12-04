function gust = isimip_windfield_holland(tc_track,centroids,~,silent_mode,~)
% TC windfield calculation
% NAME:
%   isimip_windfield_holland
% PURPOSE:
%   function to determine advanced winfdield
%
%   given a TC track (lat/lon,CentralPressure,MaxSustainedWind), calculate
%   the wind field at locations (=centroids)
%
%   See also the original (not fully functional) Matlab version of
%   helper_advanced_windfield_global.py, i.e.
%   helper_advanced_windfield_global.m
%
%   so far adapted to windfield by holland1980, holland2008, holland2010 cf peduzzi2012
%   can be easily expanded to different wind-fields eg chavez2015
%
% CALLING SEQUENCE:
%   gust = isimip_windfield_holland(tc_track,centroids,~,silent_mode,~)
%   [vfull,pcen,vmax]=isimip_windfield_holland(msize,res,cgps,ngps,penv,pcen,rmax,vmax,tint,prepcen,model)
% EXAMPLE:
%   tc_track=isimip_ibtracs_read('TEST');
%   tc_track=climada_tc_equal_timestep(tc_track); % make equal timestep
%   gust=isimip_windfield_holland(tc_track)
% INPUTS:
%   tc_track: a structure with the single track information (length(tc_track)!=1)
%       see e.g. climada_tc_read_unisys_tc_track
%       tc_track.Azimuth and/or tc_track.Celerity calculated, if not existing
%       but climada_tc_equal_timestep mist have been run and
%       tc_track.MaxSustainedWind must exist on input
%   centroids: a structure with the centroids information (see e.g.
%       climada_centroids_read):
%       centroids.lat: the latitude of the centroids
%       centroids.lon: the longitude of the centroids
% OPTIONAL INPUT PARAMETERS:
%   silent_mode: default=0, if =-1, use step-by-step detailed windfield,
%       i.e. reduce wind to zero at center of the eye (not recommended for
%       probabilistic, since hit/miss issue with closest node, see variable
%       max_wind_at_bullseye in code).
% OUTPUTS:
%   gust: the windfield [m/s] at all centroids, NOT sparse for speedup
%       i.e. convert like hazard.intensity()=sparse(res.gust)...
% RESTRICTIONS:
% MODIFICATION HISTORY:
% Tobias Geiger, geiger@pik-potsdam.de, 2016, Copyright, python verison
% David N. Bresch, david.bresch@gmail.com, 20161204
%-

gust = []; % init output

%global climada_global
% for SPEEDUP, we assume init_vars to have been executed
%if ~climada_init_vars, return; end

if ~exist('tc_track' ,'var'),return; end
if ~exist('centroids','var'),centroids=[]; end
if ~exist('silent_mode','var'),silent_mode=1; end

% PARAMETERS
%
% the model to be used
model='H08';
%
rho=1.15; % data for windfield calculation
EnvironmentalPressure_default=1010; % mb
%
% TEST area (Southern tip of Florida)
TEST_lon=-82:0.25:-80;TEST_lat=25:0.25:27;

if isempty(centroids) % generate TEST centroids
    [centroids.lon,centroids.lat] = meshgrid(TEST_lon,TEST_lat); % get 2-D
    centroids.lon=reshape(centroids.lon,[1 numel(centroids.lon)]); % make 1-D
    centroids.lat=reshape(centroids.lat,[1 numel(centroids.lat)]); % make 1-D
end

% use default EnvironmentalPressure, if not provided
if ~isfield(tc_track,'EnvironmentalPressure')
    tc_track.EnvironmentalPressure=tc_track.CentralPressure*0+EnvironmentalPressure_default;
end

if ~isfield(tc_track,'RadiusMaxWind') % set with default
    tc_track.RadiusMaxWind=tc_track.CentralPressure*0;
end

% make sure that CentralPressure never exceeds EnvironmentalPressure
% (original: if pcen > penv,pcen=penv;end with pcen=tc_track.CentralPressure)
p_points=tc_track.CentralPressure>tc_track.EnvironmentalPressure;
tc_track.CentralPressure(p_points)=tc_track.EnvironmentalPressure(p_points);

% extrapolate RadiusMaxWind from pressure if not given
% (original: see extra_rmax)
% cubic fit , ibtracs 1980-2013 (r2=0.22)
% old: rmax_from_p= -45201.105207 + 146.043463*pcen - 0.157263*pcen**2 + 0.000056*pcen**3 + 0.097515*lat + 0.016056*lon
% rmax extrapolated from historic distribution
% rmax-thresholds in nm
rmax1=15;rmax2=25;rmax3=50;
p1=950;p2=980;p3=1020;
tc_track.RadiusMaxWind(tc_track.CentralPressure<=p1)=rmax1;
p_points=tc_track.CentralPressure>p1 & tc_track.CentralPressure<=p2;
tc_track.RadiusMaxWind(p_points)=(tc_track.CentralPressure(p_points)-p1)*(rmax2-rmax1)/(p2-p1) + rmax1;
p_points=tc_track.CentralPressure>p2;
tc_track.RadiusMaxWind(p_points)=(tc_track.CentralPressure(p_points)-p2)*(rmax3-rmax2)/(p3-p2) + rmax2;
tc_track.RadiusMaxWind=tc_track.RadiusMaxWind*1.852; % nautical mile to km

n_nodes=length(tc_track.lon);

gust=centroids.lon*0; % init output

%node_i=32; % a point right within the TEST centroids

for node_i=2:n_nodes % process each node (later speedup potential)
    
    % calculate distance to all centroids
    r_arr=isimip_haversine(centroids.lon,centroids.lat,tc_track.lon(node_i),tc_track.lat(node_i));
    
    % calculate translation speed of track
    dist=isimip_haversine(tc_track.lon(node_i-1),tc_track.lat(node_i-1),tc_track.lon(node_i),tc_track.lat(node_i));
    dist= dist/1.852; % dist to nautical miles
    % hours = hours between track coordiantes
    vtrans=dist/tc_track.TimeStep(node_i); % nautical miles/hour
    if vtrans>30,vtrans=30;end % limit to 30 nmph
    
    % convert variable names to the ones used in isimip
    vmax    =tc_track.MaxSustainedWind(node_i);
    rmax    =tc_track.RadiusMaxWind(node_i);
    pcen    =tc_track.CentralPressure(node_i);
    prepcen =tc_track.CentralPressure(node_i-1);
    penv    =tc_track.EnvironmentalPressure(node_i);
    tint    =tc_track.TimeStep(node_i);
    yc      =tc_track.lat(node_i);
    
    % adjust pressure at previous track point
    if prepcen < 850,prepcen=pcen;end
    
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
    % avoid negative wind speed, set to zero
    vwind(vwind < 0)=0;
    
    % calculate vtrans wind field array assuming that effect of vtrans decreases with distance from eye
    r_arr_normed=rmax./r_arr;
    r_arr_normed(r_arr_normed>1)=1;
    vtrans_arr=vtrans*r_arr_normed;
    
    vfull=vwind;
%     % dynamic
%     vfull=dynamic_windfield(vwind,vtrans_arr,angarr,dist,cgps,ngps);

    gust=max(gust,vfull); % keep maximum instantaneous wind

end % node_i

figure('Name','isimip')
plot(tc_track.lon,tc_track.lat,'-k');
hold on
plot(tc_track.lon,tc_track.lat,'xk');
plot(centroids.lon,centroids.lat,'xr');
climada_color_plot(gust,centroids.lon,centroids.lat);
plot(tc_track.lon,tc_track.lat,'-k');
plot(tc_track.lon,tc_track.lat,'xk');
plot(centroids.lon,centroids.lat,'xr');

figure('Name','climada')
climada_gust = climada_tc_windfield(tc_track,centroids);
plot(tc_track.lon,tc_track.lat,'-k');
hold on
plot(tc_track.lon,tc_track.lat,'xk');
plot(centroids.lon,centroids.lat,'xr');
climada_color_plot(climada_gust,centroids.lon,centroids.lat);
plot(tc_track.lon,tc_track.lat,'-k');
plot(tc_track.lon,tc_track.lat,'xk');
plot(centroids.lon,centroids.lat,'xr');

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

% return
% 
% %%%%%%%%%%%%%%%%%%%
% 
% angarr=angular_array(cgps,msize);
% 
% 
% 
% 
% b=b_value(vtrans,vmax,penv,pcen,rho);
% xx=x_value(penv,pcen);
% bs=bs_value(vtrans,penv,pcen,prepcen,yc,xx,tint);
% 
% if strcmp(model,'H80')
%     vwind=stat_holland(r_arr,rmax,b ,penv,pcen,yc,rho);
% elseif strcmp(model,'H08')
%     % same as H80 but with bs instead of b, additional info required previous track point central pressure for pressure gradient
%     vwind=stat_holland(r_arr,rmax,bs,penv,pcen,yc,rho);
% elseif strcmp(model,'H08_v')
%     % same as H08 but uses max-wind speed instead of pressure
%     vwind=stat_holland_alt(r_arr,rmax,vmax,bs,yc);
% elseif strcmp(model,'H10')
%     % in holland2010 x_exp varies with distance, not implemented
%     % so far identical to holland2008 but without gradient wind
%     % note: NO gradient wind implemented for holland2010!!!
%     x_exp=0.5;
%     vwind=stat_holland2010(r_arr,rmax,penv,pcen,bs,x_exp,rho);
% end
% % avoid negative wind speed, set to zero
% vwind(vwind < 0)=0;
% 
% % calculate vtrans wind field array assuming that effect of vtrans decreases with distance from eye
% r_arr_normed=rmax./r_arr;
% r_arr_normed(r_arr_normed>1)=1;
% vtrans_arr=vtrans*r_arr_normed;
% 
% % dynamic
% vfull=dynamic_windfield(vwind,vtrans_arr,angarr,dist,cgps,ngps);
% 
% %vfull=np.round(vfull).astype(int)
% 
% end % isimip_windfield_holland
% 
% % follow helper functions
% 
% function km=haversine(lon1,lat1,lon2,lat2)
% % Calculate the great circle distance between two points
% % on the earth (specified in decimal degrees)
% 
% % convert decimal degrees to radians
% lon1=lon1/180*pi;
% lon2=lon2/180*pi;
% lat1=lat1/180*pi;
% lat2=lat2/180*pi;
% 
% % haversine formula
% dlon = lon2 - lon1;
% dlat = lat2 - lat1;
% a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
% c = 2 * asin(sqrt(a));
% 
% % 6367 km is the radius of the Earth
% km = 6367 * c;
% end % haversine
% 
% function [xf,yf]=final_coord(i,j,xc,yc,res,temp)
% xf=(j-temp)*res+xc;
% yf=(temp-i)*res+yc;
% end % final_coord
% 
% function v_arr=stat_holland(r_arr,rmax,b,penv,pcen,lat,rho)
% % holland symmetric and static wind field according to Holland1980 and with b=bs according to Holland2008
% lat=abs(lat);
% f=2*0.0000729 * sin(lat*pi/180);
% % units are m/s
% v_arr=   sqrt(((100*b/rho*(rmax./r_arr).^b*(penv-pcen)*   exp(-(rmax./r_arr).^b)).^0.5).^2+(1000*0.5*r_arr*f) ^2) - 0.5*1000*r_arr*f;
% % v_arr=np.sqrt(((100*b/rho*(rmax/r_arr)**b*(penv-pcen)*np.exp(-(rmax/r_arr)**b))**0.5)**2+(1000*0.5*r_arr*f)**2) - 0.5*1000*r_arr*f
% v_arr(isnan(v_arr))=0;
% % translate to knots
% v_arr=v_arr/0.51444444444;
% end % stat_holland
% 
% function v_arr=stat_holland_alt(r_arr,rmax,vms,b,lat)
% % holland symmetric and static wind field according to Holland1980 and with b=bs according to Holland2008
% lat=abs(lat);
% f=2*0.0000729 * sin(lat*pi/180.);
% % units are m/s
% vms=vms*0.51444444444;
% v_arr=   sqrt(vms ^2*((rmax/r_arr) ^b * exp(1-(rmax/r_arr) ^b))+(1000*0.5*r_arr*f) ^2) - 0.5*1000*r_arr*f;
% % v_arr=np.sqrt(vms**2*((rmax/r_arr)**b * exp(1-(rmax/r_arr)**b))+(1000*0.5*r_arr*f)**2) - 0.5*1000*r_arr*f
% v_arr(isnan(v_arr))=0;
% % translate to knots
% v_arr=v_arr/0.51444444444;
% end % stat_holland_alt
% 
% function v_arr=stat_holland2010(r_arr,rmax,penv,pcen,bs,x_exp,rho)
% v_arr= ((100*bs*(penv-pcen)*(rmax/r_arr) ^bs) / (rho*exp((rmax/r_arr) ^bs))) ^x_exp;
% % v_arr= ((100*bs*(penv-pcen)*(rmax/r_arr)**bs) / (rho*exp((rmax/r_arr)**bs)))**x_exp
% 
% % translate to knots
% v_arr=v_arr/0.51444444444;
% end % stat_holland2010
% 
% 
% 
% function [vtrans,dist]=v_trans_speed(lon1,lat1,lon2,lat2,tint)
% % calculate distance to next grid point and translational cyclone speed
% dist=haversine(lon1,lat1,lon2,lat2);
% dist= dist/1.852; % dist to nautical miles
% % hours = hours between track coordiantes
% hours=tint;
% vtrans=dist/hours;
% if vtrans>30,vtrans=30;end
% end % v_trans_speed
% 
% function b=b_value(vtrans,vmax,penv,pcen,rho)
% if vmax>0
%     b=rho*   exp(1)*(0.51444444444*(vmax-vtrans)) ^2/(100*(penv-pcen));
%     % b=rho*np.exp(1)*(0.51444444444*(vmax-vtrans))**2/(100*(penv-pcen))
%     if b < 1
%         b=1.0;
%     elseif b > 2.5
%         b=2.5;
%     else
%         b=1.0;
%     end
% else
%     b=1; % ADDED dnb
% end
% end % b_value
% 
% 
% function bs=bs_value(vtrans,penv,pcen,prepcen,lat,xx,tint)
% vt_ms=vtrans*0.51444444444;
% bs=-4.4e-5*(penv-pcen) ^2 + 0.01*(penv-pcen) + 0.03*(pcen-prepcen)/tint - 0.014*   abs(lat) + 0.15*vt_ms ^xx + 1.0;
% % bs=-4.4e-5*(penv-pcen)**2 + 0.01*(penv-pcen) + 0.03*(pcen-prepcen)/tint - 0.014*np.abs(lat) + 0.15*vt_ms**xx + 1.0
% end % bs_value
% 
% function v_from_p=extra_v(pcen,lat,lon)
% % determined in knots 10min sustained wind
% %v_from_p = 1131.816+0.064*lat-0.031*lon-1.1000*pcen  # peduzzi
% v_from_p = 1142.000 + 0.056*lat - 0.032*lon - 1.111*pcen;  % ibtracs 1980-2013 (r2=0.91)
% end % extra_v
% 
% function p_from_v=extra_p(vmax,lat,lon)
% % determined in hPa
% %p_from_v = 1024.688+0.055*lat-0.028*lon-0.815*vmax # peduzzi
% p_from_v = 1024.388 + 0.047*lat - 0.029*lon - 0.818*vmax; % ibtracs 1980 -2013 (r2=0.91)
% end %  extra_p
% 
% function rmax_from_p=extra_rmax(pcen) % ,lat,lon)
% % cubic fit , ibtracs 1980-2013 (r2=0.22)
% %rmax_from_p= -45201.105207 + 146.043463*pcen - 0.157263*pcen**2 + 0.000056*pcen**3 + 0.097515*lat + 0.016056*lon
% % rmax extrapolated from historic distribution
% % rmax-thresholds in nm
% rmax1=15; rmax2=25; rmax3=50;
% p1=950; p2=980; p3=1020;
% if pcen <= p1
%     val=rmax1;
% elseif  (pcen > p1) && (pcen <= p2)
%     val= (pcen-p1)*(rmax2-rmax1)/(p2-p1) + rmax1;
% elseif (pcen > p2)
%     val= (pcen-p2)*(rmax3-rmax2)/(p3-p2) + rmax2;
% end
% rmax_from_p=val;
% end % rmax_from_p
% 
% function angarr=angular_array(cgps,msize)
% % define general array that contains angular wind-speed dependence for translational wind speed of cyclone
% % northward movement of cyclone assumed
% angarr=zeros(msize,msize);
% temp=ceil(msize/2);
% %temp=msize//2
% %xc=cgps(1);
% yc=cgps(2);
% for i=1:msize
%     for j=1:msize
%         xf=j-temp;
%         yf=temp-i;
%         coss = xf/sqrt(yf^2+xf^2);
%         if yc >0 % northern hemisphere
%             if yf >= 0
%                 angarr(i,j)=acos(coss)*180/pi;
%             else
%                 angarr(i,j)=360 - acos(coss)*180/pi;
%             end
%         else % southern hemisphere
%             if yf >= 0
%                 angarr(i,j)=acos(coss)*180/pi;
%             else
%                 angarr(i,j)=360 - acos(coss)*180/pi;
%             end
%         end
%     end % j
% end %i
% 
% end % angular_array
% 
% function vfull=dynamic_windfield(vwind,vtrans,angarr,dist,cgps,ngps)
% % calculate angle of direction using vtrans and next track point
% % determine direction, northern hemisphere deviation in degree from northward movement
% 
% % center coordinates from which to measure distance
% xc=cgps(1);
% yc=cgps(2);
% % next track point
% xn=ngps(1);
% yn=ngps(2);
% 
% londiff=xn-xc;
% latdiff=yn-yc;  % if latdiff ... add or subract 90 degrees, test all possible configurations
% londist=haversine(xc,yc,xn,yc); % distance along same latitude
% londist= londist/1.852;
% %latdist=haversine(xc,yc,xc,yn);  % distance along same longitude
% %latdist= latdist/1.852;
% % londiff for values abs(londiff) > 300 happen at point where 180 turns to -179 etc
% if londiff < 0
%     if londiff < -300
%         pref=1;
%     else
%         pref=-1;
%     end
% else
%     if londiff > 300
%         pref=-1;
%     else
%         pref=1;
%     end
% end
% if dist>0
%     if latdiff < 0
%         direct=pref*(180 - asin(londist/dist)*180/pi);
%         %direct=pref*(180 - np.arccos(latdist/dist)*180/np.pi)
%     else
%         direct=pref*asin(londist/dist)*180/pi;
%         %direct=pref*np.arccos(latdist/dist)*180/np.pi
%     end
% else
%     fprintf('no movement of cyclone\n')
%     direct=0;
% end
% 
% angarr=angarr + direct;
% angarr=cos(angarr*pi/180.);
% angarr(isnan(angarr))=0;
% if yc < 0
%     %angarr=np.flipud(angarr)
%     %angarr=np.fliplr(angarr)
%     angarr=angarr'; % GUESS dnb
% end
% % calculate windfield plus translation speed
% vfull=   sqrt(vwind ^2 + vtrans ^2 + 2*vwind*vtrans*angarr);
% % vfull=np.sqrt(vwind**2 + vtrans**2 + 2*vwind*vtrans*angarr)
% end % dynamic_windfield
