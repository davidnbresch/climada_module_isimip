function hazard=isimip_flood_load(flood_fraction_filename,flood_depth_filename,hazard_filename,entity,check_plot)
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
%   next call: climada_EDS_calc
% CALLING SEQUENCE:
%   hazard=isimip_flood_load(flood_fraction_filename,flood_depth_filename,hazard_filename,entity,check_plot)
% EXAMPLE:
%   entity=climada_entity_load('USA_UnitedStates_Florida');
%   hazard=isimip_flood_load('fldfrc_max.nc','flddph_max.nc','auto',entity);
% INPUTS:
%   flood_fraction_filename: filename of the .nc file with the flood
%       fraction footprints, default folder is ..climada_data/isimip
%       > promted for if not given
%       fraction is in the range 0..1
%   flood_depth_filename: filename of the .nc file with the flood
%       depth footprints, default folder is ..climada_data/isimip
%       > promted for if not given
%       depth in units of meters [m]
%   hazard_filename: the filename (with or without path) the generated
%       hazard set is stored to. If='auto', the name is autmatically
%       generated, by appendign _FL to the entity name (still stored into
%       ../hazards folder)
%   entity: an entity struct to interpolate the flood footprints to, see
%       climada_entity_load and climada_entity_read for a description
% OPTIONAL INPUT PARAMETERS:
%   check_plot: whether show a check plot (=1, default), or not (=0)
%       Note that plotting might often take longer than the full
%       conversion...
%       if =2, also show the centroids as red dots
%       if =3, also show the original data grid as blue dots (might take time...)
% OUTPUTS:
%   hazard: a climada hazard structure, see manual
%       in addition to the standard hazard.intensity, this hazard also
%       contains hazard.fraction, the flooded fraction of the area the
%       centroid represents. See special case in climada_EDS_calc, too.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160929, initial
% David N. Bresch, david.bresch@gmail.com, 20160930, generalized, hazard.fraction added
%-

hazard=[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('flood_fraction_filename','var'),flood_fraction_filename='';end
if ~exist('flood_depth_filename','var'),   flood_depth_filename=   '';end
if ~exist('hazard_filename','var'),        hazard_filename=        '';end
if ~exist('entity','var'),                 entity=                 '';end
if ~exist('check_plot','var'),             check_plot=              1;end


% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% PARAMETERS
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end
%
sparse_density=.01; % density of hazard.intensity (sparse, guess to allocate)
%
% interpolation method, see help interpn
interpn_method='linear'; % default 'linear', also: 'nearest','spline','cubic'

% prompt for flood_fraction_filename if not given
if isempty(flood_fraction_filename) % local GUI
    flood_fraction_filename=[isimip_data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(flood_fraction_filename, 'Select flood fraction file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        flood_fraction_filename=fullfile(pathname,filename);
    end
end

% prompt for flood_depth_filename if not given
if isempty(flood_depth_filename) % local GUI
    flood_depth_filename=[isimip_data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(flood_depth_filename, 'Select flood depth file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        flood_depth_filename=fullfile(pathname,filename);
    end
end

% prompt for hazard_filename if not given
if isempty(hazard_filename) % local GUI
    hazard_filename=[climada_global.hazards_dir filesep '*.mat'];
    [filename, pathname] = uiputfile(hazard_filename, 'Save hazard as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_filename=fullfile(pathname,filename);
    end
end

% flood_fraction_filename: complete path, if missing
[fP,fN,fE]=fileparts(flood_fraction_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_fraction_filename=[fP filesep fN fE];
if ~exist(flood_fraction_filename,'file')
    fprintf('Error: %s not found\n',flood_fraction_filename);
end

% flood_depth_filename: complete path, if missing
[fP,fN,fE]=fileparts(flood_depth_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_depth_filename=[fP filesep fN fE];
if ~exist(flood_depth_filename,'file')
    fprintf('Error: %s not found\n',flood_depth_filename);
end

if strcmpi(hazard_filename,'auto') % assign automatically
    [~,fN]=fileparts(entity.assets.filename);
    hazard_filename=strrep(fN,'_entity','');
    hazard_filename=[climada_global.hazards_dir filesep hazard_filename '_FL.mat'];
end

% hazard_filename: complete path, if missing
[fP,fN,fE]=fileparts(hazard_filename);
if isempty(fP),fP=climada_global.hazards_dir;end
if isempty(fE),fE='.mat';end
hazard_filename=[fP filesep fN fE];

% load entity to obtsin centroids
entity=climada_entity_load(entity);

fprintf('reading lon, lat and time from %s and \n  %s ...',flood_fraction_filename,flood_depth_filename);
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon      = ncread(flood_fraction_filename,'lon');
nc.lat      = ncread(flood_fraction_filename,'lat');
nc.time     = ncread(flood_fraction_filename,'time');
fprintf(' done\n');

% find the bounding box around the assets
lonmin=min(entity.assets.lon);lonmax=max(entity.assets.lon);
latmin=min(entity.assets.lat);latmax=max(entity.assets.lat);

% figure grid spacing
dlon=abs(max(diff(nc.lon)));dlat=abs(max(diff(nc.lat)));

% figure the index range to cover a bot more than the box around the assets
[~,lon_index_min]=min(abs(nc.lon-(lonmin-2*dlon)));
[~,lon_index_max]=min(abs(nc.lon-(lonmax+2*dlon)));
[~,lat_index_min]=min(abs(nc.lat-(latmin-2*dlat)));
[~,lat_index_max]=min(abs(nc.lat-(latmax+2*dlat)));

fprintf('reading lon/lat %2.2f .. %2.2f/%2.2f .. %2.2f\n',nc.lon(lon_index_min),nc.lon(lon_index_max),nc.lat(lat_index_min),nc.lat(lat_index_max));

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
n_events   =length(nc.time); % number of events
fprintf('generating FL hazard set for %i events at %i centroids\n',n_events,n_centroids);
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
hazard.yyyy=(1:n_events)+hazard.reference_year-nc.time(end)-1;
hazard.mm=ones(1,n_events);
hazard.dd=ones(1,n_events);
hazard.datenum=datenum(hazard.yyyy,hazard.mm,hazard.dd);
hazard.scenario='no clima;te change';
hazard.name=cellstr(num2str(hazard.event_ID'))';
hazard.frequency=ones(1,n_events)/n_events;
hazard.filename=hazard_filename;
hazard.date=datestr(now);
hazard.comment=sprintf('FL hazard event set, generated %s',hazard.date);
% and the sparse arrays:
hazard.intensity=spalloc(n_events,n_centroids,ceil(n_events*n_centroids*sparse_density));
hazard.fraction =spalloc(n_events,n_centroids,ceil(n_events*n_centroids*sparse_density));

% loop over each event, get depth (anf fraction) and store to centroids
% template for-loop with waitbar or progress to stdout
t0       = clock;
mod_step = 1; % first time estimate after 10 events, then every 100
format_str='%s';

for event_i=1:n_events
    
    % get single timestep and reduced 'tile'
    nc.fraction = ncread(flood_fraction_filename,'var',[lon_index_min lat_index_min event_i],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    nc.depth    = ncread(flood_depth_filename,   'var',[lon_index_min lat_index_min event_i],[index_dlon+1 index_dlat+1 1]); % lon, lat, time
    
    % interpolate to centroids
    depth_tile   =interpn(tile_lon,tile_lat,nc.depth,hazard.lon,hazard.lat,interpn_method);
    depth_tile(isnan(depth_tile))      =0; % replace NaN with zeros
    if isfield(hazard,'isgridpoint'),depth_tile(hazard.isgridpoint)   =0;end % regular grid 'around' assets
    hazard.intensity(event_i,:)=depth_tile;
    
    fraction_tile=interpn(tile_lon,tile_lat,nc.fraction,hazard.lon,hazard.lat,interpn_method);
    fraction_tile(isnan(fraction_tile))=0; % replace NaN with zeros
    if isfield(hazard,'isgridpoint'),fraction_tile(hazard.isgridpoint)=0;end % regular grid 'around' assets
    hazard.fraction(event_i,:) =fraction_tile;
    
    % the progress management
    if mod(event_i,mod_step)==0
        mod_step          = 10;
        t_elapsed_event   = etime(clock,t0)/event_i;
        events_remaining  = n_events-event_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('converting ..., est. %3.0f sec left (%i/%i events)',t_projected_sec,   event_i,n_events);
        else
            msgstr = sprintf('converting ..., est. %3.1f min left (%i/%i events)',t_projected_sec/60,event_i,n_events);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end
fprintf(format_str,''); % move carriage to begin of line

hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);

fprintf('saving hazard as %s\n',hazard_filename)
save(hazard_filename,'hazard');

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
        [X,Y] = meshgrid(tile_lon,tile_lat);hold on;plot(X,Y,'.r','MarkerSize',1);
        fprintf(' done\n');
        legend({'hazard intensity','centroids','data grid'});
    end % check_plot
end % plot max hazard at each centroid

end % isimip_flood_load