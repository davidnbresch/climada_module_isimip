function hazard=isimip_flood_load(flood_filename,hazard_filename,entity,check_plot,isimip_simround,years_range,interpn_method,silent_mode,subtract_matsiro)

% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_flood_load
% PURPOSE:
%   load .mat file with flood footprints and construct a climada flood
%   hazard event set
%
%   flood data reference: sven.willner@pik-potsdam.de
%
%   Sepcial: since flood footrpints are 'tiny', i.e. should not be
%   displayed using interpolation between far distand grid points, this
%   code does check for entity.assets.isgridpoint and only assigns flood
%   depth and fraction values to 'real' asset centroids, not to the grid
%   'around' (see code)  
%
%   See isimip_flood_import to just the country ISO3 code to define the
%   geographic area of the hazard set. 
%
%   Octave: please install the netCDF package first: 
%   pkg install -forge netcdf -local -auto
%   or: pkg install -verbose -forge -auto netcdf
%
%   next call: climada_EDS_calc
% CALLING SEQUENCE:
%   hazard=isimip_flood_load(flood_filename,hazard_filename,entity,check_plot)
% EXAMPLE:
%   entity=climada_entity_load('USA_UnitedStates_Florida');
%   hazard=isimip_flood_load('USA_UnitedStates_Florida_FL.nc','auto',entity);
%   hazard=isimip_flood_load('global_FL.nc','auto',entity);
% INPUTS:
%   flood_filename: filename of the .nc file with the flood
%       footprints, default folder is ..climada_data/isimip/FL
%       > promted for if not given
%       fraction (variable name 'fldfrc') is in the range 0..1
%       depth (variable name 'flddph') in units of meters [m]
%       there should be one event per year (i.e., yearly maxima)
%   hazard_filename: the filename (with or without path) the generated
%       hazard set is stored to. If='auto', the name is automatically
%       generated using function 'isimip_get_flood_hazard_filename'
%   entity: an entity struct to interpolate the flood footprints to, see
%       climada_entity_load and climada_entity_read for a description
% OPTIONAL INPUT PARAMETERS:
%   check_plot: whether show a check plot (=1, default), or not (=0)
%       Note that plotting might often take longer than the full
%       conversion...
%       if =2, also show the centroids as red dots
%       if =3, also show the original data grid as blue dots (might take time...)
%   isimip_simround: the sub-directory within the isimip folder within
%       climada_data/isimip where the raw isimip data (NetCDF file
%       flood_filename) is located. This does not affect where the hazard is
%       saved. If not specified, the file is assumed to be located directly
%       within the climada_data/isimip folder.
%   years_range: vector of length 2 containing the first and the last year
%       to be loaded from the netcdf file. If empty or [0 0], loads all data.
%   interpn_method: interpolation method to map hazard onto the entity
%       structure (i.e., centroids). The default value is 'linear' and
%       should be used if the hazard it at a lower resolution than the
%       entity, however it should be set to 'nearest' when using centroids
%       on the same grid (e.g., from a file at a as0150 resolution such as
%       'gdp_1980-2100_SSP2_0150as_remapnn_yearly.nc'. Possible values are:
%       'linear', 'nearest','spline','cubic'.
%   silent_mode: if 1, no output is given (expect errors etc), if 0 many
%       output given
%   subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%       fraction from the data. Default =0.
% OUTPUTS:
%   hazard: a climada hazard structure, see manual
%       in addition to the standard hazard.intensity, this hazard also
%       contains hazard.fraction, the flooded fraction of the area the
%       centroid represents. See special case in climada_EDS_calc, too.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160929, initial
% David N. Bresch, david.bresch@gmail.com, 20160930, generalized, hazard.fraction added
% Sven Willner, sven.willner@pik-potsdam.de, 20160930, reduced to one flood file
% David N. Bresch, david.bresch@gmail.com, 20161002, small fix to print netCDF filename in stdout
% David N. Bresch, david.bresch@gmail.com, 20170705, climada_global.save_file_version
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20171130, add optional
%   argument 'isimip_data_subdir', and improved treatment of the time axis
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180118, fixes to use
% ISIMIP-2a (previous version kept as isimip_flood_load_orig.m
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180518, new function
%   argument 'subtract_matsiro'
% David N. Bresch, david.bresch@gmail.com, 20190211, subtract_matsiro default set
% David N. Bresch, david.bresch@gmail.com, 20190309, hint to isimip_flood_import
%-

hazard=[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('flood_filename','var'),         flood_filename=         '';end
if ~exist('hazard_filename','var'),        hazard_filename=        '';end
if ~exist('entity','var'),                 entity=                 '';end
if ~exist('check_plot','var'),             check_plot=              1;end
if ~exist('years_range','var'),            years_range=         [0 0];end
if ~exist('interpn_method','var')          interpn_method=   'linear';end
if ~exist('isimip_simround','var')         isimip_simround=   '';end
if ~exist('silent_mode','var')             silent_mode=0;end
if isequal(isimip_simround, '')
    isimip_data_dir = [climada_global.data_dir filesep 'isimip'];
else
    isimip_data_dir = [climada_global.data_dir filesep 'isimip' filesep isimip_simround];
end
if ~exist('subtract_matsiro','var')        subtract_matsiro=0;end

% check validity of arguments
if ~isequal(size(years_range), [1 2])
    warning('** error ** years_range should be of size [1 2] *****')
    return
end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% PARAMETERS
%
%
sparse_density=.01; % density of hazard.intensity (sparse, guess to allocate)

% check that flood_filename is given
if isequal(flood_filename,'') % local GUI
    fprintf('Error: flood_filename is not given')
    return; % cancel
end

% check that hazard_filename is given
if isequal(hazard_filename,'') % local GUI
    fprintf('Error: hazard_filename is not given')
    return; % cancel
end

% flood_filename: complete path, if missing
[fP,fN,fE]=fileparts(flood_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_filename=[fP filesep fN fE];
if ~exist(flood_filename,'file')
    fprintf('Error: %s not found\n',flood_filename);
end

if strcmpi(hazard_filename,'auto') % assign automatically
    hazard_filename=isimip_get_flood_hazard_filename(flood_filename,entity,isimip_simround,years_range);
%     [~,fN]=fileparts(entity.assets.filename);
%     [~,fN2]=fileparts(flood_filename);
%     hazard_filename=strrep(fN,'_entity','');
%     if ~isequal(years_range, [0 0])
%         hazard_filename=[climada_global.hazards_dir filesep fN2 '_' fN '_FL.mat'];
%     else
%         hazard_filename=[climada_global.hazards_dir filesep fN2 '_' fN '_' years_range(1) '-' years_range(2) '_FL.mat'];
%     end
end
% 
% % hazard_filename: complete path, if missing
% [fP,fN,fE]=fileparts(hazard_filename);
% if isempty(fP),fP=climada_global.hazards_dir;end
% if isempty(fE),fE='.mat';end
% hazard_filename=[fP filesep fN fE];

% load entity to obtsin centroids
entity=climada_entity_load(entity);


if ~silent_mode,fprintf('reading lon, lat and time from %s ...',flood_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon      = ncread(flood_filename,'lon');
nc.lat      = ncread(flood_filename,'lat');
nc.time     = ncread(flood_filename,'time');
nc.time_units = ncreadatt(flood_filename,'time','units');

% Get info about NetCDF time axis (used later)
n_events_orig   =length(nc.time); % number of events
nc.time_units = strsplit(nc.time_units, ' ');
nc.time_orig = cell2mat(strcat(nc.time_units(3), {' '}, nc.time_units(4)));
nc.time_units = cell2mat(nc.time_units(1));
if ~silent_mode,fprintf(' done\n');end

% find the bounding box around the assets
lonmin=min(entity.assets.lon);lonmax=max(entity.assets.lon);
latmin=min(entity.assets.lat);latmax=max(entity.assets.lat);

% figure grid spacing
dlon=max(abs(diff(nc.lon)));dlat=max(abs(diff(nc.lat)));

% figure the index range to cover a bit more than the box around the assets
[~,lon_index_min]=min(abs(nc.lon-(lonmin-2*dlon)));
[~,lon_index_max]=min(abs(nc.lon-(lonmax+2*dlon)));
[~,lat_index_min]=min(abs(nc.lat-(latmin-2*dlat)));
[~,lat_index_max]=min(abs(nc.lat-(latmax+2*dlat)));

if ~silent_mode
    fprintf('reading lon/lat %2.2f .. %2.2f/%2.2f .. %2.2f\n',nc.lon(lon_index_min),nc.lon(lon_index_max),nc.lat(lat_index_min),nc.lat(lat_index_max));
    fprintf('       (indices %d .. %d/%d .. %d)\n',lon_index_min,lon_index_max,lat_index_min,lat_index_max);
end

lon_index_temp=lon_index_min; % can be ordered descending
lon_index_min=min(lon_index_min,lon_index_max);
lon_index_max=max(lon_index_temp,lon_index_max);

lat_index_temp=lat_index_min; % can be ordered descending
lat_index_min=min(lat_index_min,lat_index_max);
lat_index_max=max(lat_index_temp,lat_index_max);

% fiddle with tile (sub-section of global file)
index_dlon=lon_index_max-lon_index_min;
index_dlat=lat_index_max-lat_index_min;

tile_lon=nc.lon(lon_index_min:lon_index_max);
tile_lat=nc.lat(lat_index_min:lat_index_max);

% populate the fields in hazard
n_centroids=length(entity.assets.lon); % number of centroids
% create dates for events
if strcmp(nc.time_units, 'years')
    hazard.mm=ones(1,n_events_orig);
    hazard.dd=ones(1,n_events_orig);
    hazard.yyyy=str2double(nc.time_orig(1:4))+nc.time';
    hazard.datenum=datenum(hazard.yyyy,hazard.mm,hazard.dd);
else
    if strcmp(nc.time_units, 'seconds')
        time_ratio = 3600*24;
    elseif strcmp(nc.time_units, 'minutes')
        time_ratio = 60*24;
    elseif strcmp(nc.time_units, 'hours')
        time_ratio = 24;
    else
        warning('** the time units in the NetCDF file are unexpected - assuming one every 366 days (roughly one year) ending in the reference year **')
        time_ratio=1/366;
        nc.time=(1:n_events_orig)-1;
    end
    nc.time_datenum = double(nc.time/time_ratio+datenum(nc.time_orig));
    hazard.datenum=nc.time_datenum;
end
% subset of events within years_range
hazard.yyyy=str2num(datestr(hazard.datenum, 'yyyy'));
if ~isequal(years_range, [0 0])
    event_keep=(hazard.yyyy>=years_range(1) & hazard.yyyy<=years_range(2));
    event_keep_which=find(event_keep);
    if ~silent_mode,fprintf('keeping subset of years (%i out of %i)',sum(event_keep),length(hazard.yyyy));end
    hazard.datenum=hazard.datenum(event_keep);
else
    event_keep=ones(1,length(hazard.yyyy));
    event_keep_which=find(event_keep);
end
hazard.yyyy=datestr(hazard.datenum, 'yyyy');
hazard.mm=datestr(hazard.datenum, 'mm');
hazard.dd=datestr(hazard.datenum, 'dd');
n_events=length(hazard.yyyy);
if ~silent_mode,fprintf('generating FL hazard set for %i events at %i centroids\n',n_events,n_centroids);end
hazard.peril_ID='FL';
hazard.units='m';
hazard.reference_year=climada_global.present_reference_year;
hazard.lon=entity.assets.lon;
hazard.lat=entity.assets.lat;
if isfield(entity.assets,'isgridpoint'),hazard.isgridpoint=entity.assets.isgridpoint;end
hazard.centroid_ID=1:n_centroids;
hazard.orig_years=n_events;
hazard.orig_event_count=n_events;
hazard.event_count=n_events;
hazard.event_ID=1:n_events;
hazard.orig_event_flag=ones(1,n_events);

% 
% hazard.yyyy=(1:n_events)+hazard.reference_year-nc.time(end)-1;
% hazard.yyyy
% hazard.mm=ones(1,n_events);
% hazard.dd=ones(1,n_events);
% hazard.datenum=datenum(hazard.yyyy,hazard.mm,hazard.dd);
hazard.scenario='no climate change';
hazard.name=cellstr(num2str(hazard.event_ID'))';
hazard.frequency=ones(1,n_events)/n_events;
hazard.filename=hazard_filename;
hazard.date=datestr(now);
hazard.comment=sprintf('FL hazard event set, generated %s',hazard.date);
% and the sparse arrays:
hazard.intensity=spalloc(n_events,n_centroids,ceil(n_events*n_centroids*sparse_density));
hazard.fraction =spalloc(n_events,n_centroids,ceil(n_events*n_centroids*sparse_density));

% loop over each event, get depth (and fraction) and store to centroids
% template for-loop with waitbar or progress to stdout
t0       = clock;
mod_step = 1; % first time estimate after 10 events, then every 100
format_str='%s';

% load 2-yr flood fraction from matsiro if needed
if subtract_matsiro
    matsiro_file = [climada_global.data_dir filesep 'isimip' filesep 'matsiro_2yr_data' filesep 'fldfrc24_2.nc'];
    matsiro_fraction = ncread(matsiro_file,'fldfrc',[lon_index_min lat_index_min 1],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    matsiro_fraction   =interpn(tile_lon,tile_lat,matsiro_fraction,hazard.lon,hazard.lat,interpn_method);
    matsiro_fraction(isnan(matsiro_fraction))      =0; % replace NaN with zeros
end

for event_i=1:n_events
    
    event_i_nc=event_keep_which(event_i);
    % get single timestep and reduced 'tile'
    nc.depth    = ncread(flood_filename,'flddph',[lon_index_min lat_index_min event_i_nc],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    try
        % try to load the flood fraction from the same file first
        nc.fraction = ncread(flood_filename,'fldfrc',[lon_index_min lat_index_min event_i_nc],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    catch exception
        % it it doesn't work, try the corresponding file with 'fldfrc'
        % instead of 'flddph' in the file name.
        nc.fraction = ncread(strrep(flood_filename,'flddph','fldfrc'),'fldfrc',[lon_index_min lat_index_min event_i_nc],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    end
    % interpolate to centroids
    depth_tile   =interpn(tile_lon,tile_lat,nc.depth,hazard.lon,hazard.lat,interpn_method);
    depth_tile(isnan(depth_tile))      =0; % replace NaN with zeros
    if isfield(hazard,'isgridpoint'),depth_tile(hazard.isgridpoint)   =0;end % regular grid 'around' assets
    hazard.intensity(event_i,:)=depth_tile;
    
    fraction_tile=interpn(tile_lon,tile_lat,nc.fraction,hazard.lon,hazard.lat,interpn_method);
    fraction_tile(isnan(fraction_tile))=0; % replace NaN with zeros
    % subtract 2-yr matsiro flood fraction if needed
    if subtract_matsiro
        fraction_tile = fraction_tile - matsiro_fraction;
        fraction_tile(fraction_tile<0)=0;
    end
    if isfield(hazard,'isgridpoint'),fraction_tile(hazard.isgridpoint)=0;end % regular grid 'around' assets
    hazard.fraction(event_i,:) =fraction_tile;
    
    % the progress management
    if mod(event_i,mod_step)==0
        mod_step          = 10;
        t_elapsed_event   = etime(clock,t0)/event_i;
        events_remaining  = n_events-event_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('converting ... est. %3.0f sec left (%i/%i events)',t_projected_sec,   event_i,n_events);
        else
            msgstr = sprintf('converting ... est. %3.1f min left (%i/%i events)',t_projected_sec/60,event_i,n_events);
        end
        if ~silent_mode,fprintf(format_str,msgstr);end % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end
if ~silent_mode,fprintf(format_str,'');end % move carriage to begin of line

hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);

if ~silent_mode,fprintf('saving hazard as %s\n',hazard_filename);end
save(hazard_filename,'hazard',climada_global.save_file_version);

if check_plot
    if check_plot==1
        climada_hazard_plot(hazard,0);
    else
        climada_hazard_plot(hazard,0,'',[],1); % plot the centroids
        %hold on;plot(hazard.lon,hazard.lat,'.r','MarkerSize',1); % plot the centroids
    end % check_plot
    if check_plot>2
        % also plot the original grid of the nc data file
        fprintf('adding original data grid (takes some time ...');
        [X,Y] = meshgrid(tile_lon,tile_lat);hold on;plot(X,Y,'.r','MarkerSize',0.01);
        fprintf(' done\n');
        legend({'hazard intensity','centroids','data grid'});
    end % check_plot
end % plot max hazard at each centroid

end % isimip_flood_load