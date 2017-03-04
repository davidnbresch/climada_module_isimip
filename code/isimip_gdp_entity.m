function [entity,params]=isimip_gdp_entity(ISO3,params)
% climada isimip entity population gdp
% MODULE:
%   isimip
% NAME:
%   isimip_gdp_entity
% PURPOSE:
%   create the entity based on isimip GDP (and add population data).
%
%   Reads the GDP grid (per time, usually once per year) from
%   val_filename, multiplies with the conversion_factor (per gridcell)
%   from con_filename and stores into climada
%   entity.assets.Values(i,j) for year i and centroid j. The code maps to
%   the GLOBAL centroids, such that one can generate one (global) hazard
%   set and all single/multi country entities run properly (see
%   entity.assets.centroid_index on otuput below).
%
%   If the matching population netCDF file exists (pop_filename), the variable
%   entity.assets.Population(i,j) contains the populaiton (number of people
%   per centroid) for year i and centroid j.  Note that year i matches with
%   year i of entity.assets.Values, i.e. Population is only stored for
%   years with Values.
%
%   Please be PATIENT the first time you run this code, as it generates the
%   global reference grid, where calculation of distance to coast for all
%   land points does take (substantial) time. See PARAMETERS to switch
%   distance_to_coast off. In the first run, it also saves a centroids
%   file, in case more than 4 mio centroids, it also saves a reduced
%   verison (see code for details).
%
%   single/multi country mode (e.g. ISO3='DEU', ISO3={'DEU','FRA'}):
%    Checks the country ID(s) (NatID) on NatID_filename and takes all gridcells
%    within the requested country(s). If the resolution of the NatID does not
%    match the population data or if there is no NatID_filename provided,
%    the code uses the climada country shape files to select the gridcells
%    within the country(s) (the code notifies to stdout).
%
%   all country mode (ISO3='ALL_IN_ONE'):
%    same as for single country, but create one global entity with
%    additional fields as described in OUTPUTS below.
%
%   Octave: please install the netCDF package first:
%    pkg install -forge netcdf -local -auto
%    or: pkg install -verbose -forge -auto netcdf
%
%   next call: isimip...
% CALLING SEQUENCE:
%   entity=isimip_gdp_entity(ISO3,params)
% EXAMPLE:
%   entity=isimip_gdp_entity('DEU') % single country entity
%   entity=isimip_gdp_entity({'DEU','FRA','ITA'}) % multi country entity
%   entity=isimip_gdp_entity('ALL_IN_ONE') % one global entity
%   params=isimip_gdp_entity('params') % get default parameters
%   % to check population for year 2030:
%   params.plot_population=1;params.year=2030;climada_entity_plot(entity,2,params)
%
%   % store the hazard also into the entity for special isimip use
%   params.hazard_file='GLB_0360as_TC_hist';
%   entity=isimip_gdp_entity('BRB',params)
% INPUTS:
% INPUTS:
%   ISO3: the ISO3 country code, e.g. 'DEU' or 'FRA' or {'DEU','FRA'} to
%       combine two countries in one entity (the entity still knows which
%       centroids belong to which country, see entity.assets.NatID)
%       if ='all': process all countries (be careful, test a few first),
%           create one single entity for each country (takes time and memory,
%           hence see also nexst option).
%       if ='ALL_IN_ONE', one entity with all countries is created
%           This is a bit a special case, as the entity is mainly to be
%           used within the ismip context.
%       > promted for (to select from a list, also multiple) if not given
%       if ='params', just return default parameters, in entity, i.e. the
%           first output, already.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    val_filename: filename of the .nc file with the (annual) GDP values
%       ='0360as': use 0.1 degree default file (default, set in PARAMETERS)
%       ='0150as': use 2.5 minutes default file (set in PARAMETERS)
%       If only a filename (without path) is passed, isimip_data_dir (set
%       in PARAMETERS) is prepended
%    pop_filename: filename of the .nc file with the population values
%       if not existing, population is not written into the entity. Same
%       conventions as for val_filename apply
%    pop2_filename: filename of the .nc file with the second set of
%       population values. if not existing, population is not written into
%       the entity. Same conventions as for pop_filename apply.
%    con_filename: filename of the .nc file with the conversion
%       factor from either GDP or population to assets (per gridpoint)
%       if val_filename is one of the short names, con_filename is
%       set to the corresponding file, if empty on input. Default set in
%       PARAMETERS. If only a filename (without path) is passed,
%       isimip_data_dir (set in PARAMETERS) is prepended.
%    NatID_filename: filename of the .nc file with the national grid IDs (to
%       assign grid cells to countries). Default set in PARAMETERS
%       if val_filename is one of the short names, NatID_filename is
%       set to the corresponding file, if empty on input
%       To avoid using NatID and just use the grid cells within a country
%       shape, set NatID_filename='ignore', but if ISO3='ALL_IN_ONE', NatID
%       is needed, hence ='ignore' does not work.
%       If only a filename (without path) is passed, isimip_data_dir (set
%       in PARAMETERS) is prepended .
%    check_plot: whether show a check plot (=1, default), or not (=0)
%    distance_to_coast: whether we calculate distant to coast (in km) of
%       call centroids (default=1), which speeds up later climada calculations
%       substantially (as coastal hazards need not to be evaluated at
%       inner-continental points). Set =0 in special cases, as initial
%       calculation might easily take 1h (we're talking about millions of
%       centroids... but since climada_distance2coast_km listens to
%       climada_global.parfor, set climada_global.parfor=1 for speedup).
%    hazard_file: if the name of a valid hazard event set, join the hazard
%       intensity as field entity.hazard.intensity(event_i,centroids) and
%       the field entity.hazard.yyyy(event_i) to contain the year of the
%       event_i. Since we use climada_hazard_load, hazard_file does neither
%       need to contain the path nor the .mat extension, e.g.
%       ='GLB_0360as_TC_hist' is ok.
%    hazard_match: if=1, match assets (population) and (historic) hazard, 
%       i.e. only keep data for years we have all three. Default=1 (keeps 
%       size manageable)
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
%   PLUS the fields
%       entity.assets.Values(n_times,n_centroids) with the assets as in the
%           variable as on the gdp file for each timestep at each centroid.
%           These values are scaled by conversion_factor, see description.
%       entity.assets.Values_yyyy: the year the Values(i,:) are valid for
%       entity.assets.Population(n_times,n_centroids) with the population
%           as in the pop file for each timestep at each centroid.
%       entity.assets.Population2(n_times,n_centroids) with the population
%           as on the second  pop file for each timestep at each centroid.
%       entity.assets.NatID: the isimip country number for each centroid
%       entity.assets.ISO3_list: the list linking isimip country numbers
%           in ISO3_list(:,2) with ISO names in ISO3_list(:,1)
%       entity.assets.centroid_index: is already pointing to the
%           correct centroids in the corresponding global centroids file
%           (as retunred in params.centroid_file), but the hazard is not
%           set yet (instead of re-encoding, just define
%           entity.assets.hazard yourself).
%       entity.hazard: the hazard set just for the centroids as in assets
%
%   for direct isimip usage, you might mainly consider the fields
%       note that we have assets at centroids for each year, hazard events within years
%       entity.assets.lon(centroid_i): the longitude of each centroid
%       entity.assets.lat(centroid_i): the latitude of each centroid
%       entity.assets.Values(year_i,centroid_i): asset value for year_i at centroid_i
%       entity.assets.Values_yyyy(year_i): the years
%       entity.assets.Population(year_i,centroid_i): population for year_i at centroid_i
%       entity.hazard.intensity(event_i,centroid_i): the hazard intensity for event_i at centroid_i
%       entity.hazard.yyyy(event_i): the year of event_i, e.g.
%           pos=find(entity.hazard.yyyy==entity.assets.Values_yyyy(1)) gives all events in first year of assets
%       entity.hazard.name{event_i}: the ibtracsID (such as 1950058S20114)
%       entity.hazard.ID_no(event_i): ibtracsID as number, such as 1971275N10176 -> 1971275.10176
%
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161017, initial
% David N. Bresch, david.bresch@gmail.com, 20161120, ALL_IN_ONE added
% David N. Bresch, david.bresch@gmail.com, 20161121, list of ISO3 allowed
% David N. Bresch, david.bresch@gmail.com, 20161123, checked for 2.5min
% David N. Bresch, david.bresch@gmail.com, 20170205, renamed to isimip_gdp_entity, overhaul, con_factor, params
% David N. Bresch, david.bresch@gmail.com, 20170209, 0150as files updated
% David N. Bresch, david.bresch@gmail.com, 20170215, ISO3 sorted, if a list
% David N. Bresch, david.bresch@gmail.com, 20170217, pop_filename, hazard_file
% David N. Bresch, david.bresch@gmail.com, 20170224, isimip_NatID_RegID instead of isimip_ISO3_list
% David N. Bresch, david.bresch@gmail.com, 20170225, filename does contain resolution, such as 0360as
% David N. Bresch, david.bresch@gmail.com, 201702304, pop2_* added
%-

entity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('ISO3','var'),   ISO3   = '';end
if ~exist('params','var'), params = struct;end

% check for some parameter fields we need
if ~isfield(params,'val_filename'),       params.val_filename='';end
if ~isfield(params,'con_filename'),       params.con_filename='';end
if ~isfield(params,'pop_filename'),       params.pop_filename='';end
if ~isfield(params,'pop2_filename'),      params.pop2_filename='';end
if ~isfield(params,'NatID_filename'),     params.NatID_filename='';end
if ~isfield(params,'check_plot'),         params.check_plot=[];end
if ~isfield(params,'distance_to_coast'),  params.distance_to_coast=[];end
if ~isfield(params,'centroids_file'),     params.centroids_file='';end % output only
if ~isfield(params,'hazard_file'),        params.hazard_file='';end
if ~isfield(params,'hazard_match'),       params.hazard_match=[];end

% PARAMETERS
%
verbose=1; % default=1, to suppress output to stdout later
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end
%
% define default filenames for population or gdp, conversion and country number (NatID)
%val_filename0360as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0360as_yearly_zip.nc'];val_variable_name='var1';
val_filename0360as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0360as_remapnn_yearly.nc'];
val_variable_name='gdp_grid';
con_filename0360as=[isimip_data_dir filesep 'GDP2Asset_converter_0360as_adv_1.nc'];
con_variable_name='conversion_factor';
pop_filename0360as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0360as_yearly_zip.nc'];
pop2_filename0360as =[isimip_data_dir filesep 'hyde_ssp2_1860-2015_0360as_yearly_zip.nc4']; % pop_variable_name='gdp_grid';
pop_variable_name='var1';
NatID_filename0360as=[isimip_data_dir filesep 'NatID_grid_0360as_adv_1.nc'];
%
%val_filename0150as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0150as_yearly_zip.nc4']; % val_variable_name='var1';
val_filename0150as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0150as_remapnn_yearly.nc']; % val_variable_name='var1';
con_filename0150as=[isimip_data_dir filesep 'GDP2Asset_converter_0150as.nc']; % con_factor
pop_filename0150as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0150as_yearly_zip.nc']; % pop_variable_name='gdp_grid';
pop2_filename0150as =[isimip_data_dir filesep 'hyde_ssp2_1860-2015_0150as_yearly_zip.nc4']; % pop_variable_name='gdp_grid';
NatID_filename0150as=[isimip_data_dir filesep 'NatID_grid_0150as.nc'];
%
centroids_file=[climada_global.data_dir filesep 'centroids' filesep 'GLB_']; % Nat_id_grid name will be appended
%
% the entity template to populate a defaiult entity
entity_template=[climada_global.entities_dir filesep 'entity_template'];
%
% admin0 shape file (for fallback selection option):
admin0_shape_file=climada_global.map_border_file;
%
% to TEST the code, i.e. reads only a few (3) timesteps
TEST_mode=0; % default=0
%
% populate default parameters in params
if isempty(params.check_plot),        params.check_plot=0;end
if isempty(params.distance_to_coast), params.distance_to_coast=1;end
if isempty(params.hazard_match),      params.hazard_match=1;end

if strcmpi(ISO3,'params'),entity=params;return;end % special case, return the full params structure

if isempty(params.val_filename),params.val_filename='0360as';end

grid_resolution_str='';
if strcmpi(params.val_filename,'0360as')
    params.val_filename  = val_filename0360as;
    params.pop_filename = pop_filename0360as;
    params.pop2_filename = pop2_filename0360as;
    if isempty(params.con_filename),params.con_filename=con_filename0360as;end
    if isempty(params.NatID_filename),params.NatID_filename=NatID_filename0360as;end
    grid_resolution_str='0360as';
elseif strcmpi(params.val_filename,'0150as')
    params.val_filename  = val_filename0150as;
    params.pop_filename = pop_filename0150as;
    params.pop2_filename = pop2_filename0150as;
    if isempty(params.con_filename),params.con_filename=con_filename0150as;end
    if isempty(params.NatID_filename),params.NatID_filename=NatID_filename0150as;end
    grid_resolution_str='0150as';
end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.val_filename,'file'),params.val_filename=[isimip_data_dir filesep params.val_filename];end
if ~exist(params.pop_filename,'file'),params.pop_filename=[isimip_data_dir filesep params.pop_filename];end
if ~exist(params.pop2_filename,'file'),params.pop2_filename=[isimip_data_dir filesep params.pop2_filename];end
if ~exist(params.con_filename,'file'),params.con_filename=[isimip_data_dir filesep params.con_filename];end
if ~exist(params.NatID_filename,'file'),params.NatID_filename=[isimip_data_dir filesep params.NatID_filename];end

% read the population or GDP data
% read only one slide (timestep) to save memory, read slab-by-slab below
% -------------------------------

% first, figure the variable name for the data
nc.info=ncinfo(params.val_filename);

% figure value units
Values_factor=1;
try
    if strcmpi(nc.info.Variables(4).Attributes(2).Value,'GDP PPP 2005 billion $ US')
        Values_factor=1e9;
    end
catch
    if verbose,fprintf('WARNING: Determining Value unit from nc failed\n');end
end % try
if verbose,fprintf('Values factor f=%g (climada asset value = f * original value from file)\n',Values_factor);end

for i=1:length(nc.info.Variables)
    if strcmpi(nc.info.Variables(i).Name,'var1'),val_variable_name='var1';end
    if strcmpi(nc.info.Variables(i).Name,'gdp_grid'),val_variable_name='gdp_grid';end
end % i

if verbose,fprintf('reading lon, lat, time and %s from %s ...',val_variable_name,params.val_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon       = ncread(params.val_filename,'lon')';
nc.lat       = ncread(params.val_filename,'lat')';
nc.time      = ncread(params.val_filename,'time')';
n_times      = length(nc.time);
if TEST_mode,n_times=min(3,n_times);nc.time=nc.time(1:n_times);end % TEST mode, read only few times
temp_data      = ncread(params.val_filename,val_variable_name,[1 1 1],[Inf Inf 1]); % only one time slab
if verbose,fprintf(' done\n');end

% make sure time is useful
time_val_yyyy=zeros(1,n_times); % init
for i=1:n_times
    time_val_yyyy(i)=str2double(datestr(nc.time(i)+datenum(1950,1,1)+1,'yyyy'));
end
if time_val_yyyy(1)<1860 || time_val_yyyy(1)>2100 % strange content of nc.time
    time_val_yyyy=1980+(1:n_times);
    fprintf('WARNING: years hard wired to %i..%i\n',time_val_yyyy(1),time_val_yyyy(end))
end

if verbose,fprintf('reading %s from %s',con_variable_name,params.con_filename);end
nc.con_factor = ncread(params.con_filename,con_variable_name);
if sum(abs(size(temp_data)-size(nc.con_factor)))>0
    fprintf('\nError: sizes of %s and %s do not match, aborted\n',val_variable_name,con_variable_name);
    return
else
    fprintf(', min/max: %2.2g/%2.2g\n',min(min(nc.con_factor)),max(max(nc.con_factor)));
end

if exist(params.pop_filename,'file')
    if verbose,fprintf('reading population from %s\n',params.pop_filename);end
    nc.time_pop=ncread(params.pop_filename,'time')';
    
    % currently hard-wired (as not properly defined in netCDF (there time just = 1..n)
    time_pop_yyyy=1860-1+(1:length(nc.time_pop));
    fprintf('WARNING: population years hard wired to %i..%i\n',time_pop_yyyy(1),time_pop_yyyy(end))
else
    time_pop_yyyy=[];
end

if exist(params.pop2_filename,'file')
    if verbose,fprintf('reading second population from %s\n',params.pop2_filename);end
    nc.time_pop2=ncread(params.pop2_filename,'time')';
    
    % currently hard-wired (as not properly defined in netCDF (there time just = 1..n)
    time_pop2_yyyy=1860-1+(1:length(nc.time_pop2));
    fprintf('WARNING: second population years hard wired to %i..%i\n',time_pop2_yyyy(1),time_pop2_yyyy(end))
else
    time_pop2_yyyy=[];
end

NatID_RegID=isimip_NatID_RegID; % get the mapping ISO3 - isimip country code

% create the grid
if verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon,gridlat] = meshgrid(nc.lon,nc.lat);
gridlon=gridlon';gridlat=gridlat'; % still as grid
vectlon=reshape(gridlon,[1 numel(gridlon)]); % as 1-D vect
vectlat=reshape(gridlat,[1 numel(gridlat)]);
if verbose,fprintf(' done\n');end

if exist(params.NatID_filename,'file')
    nc.NatIDGrid = ncread(params.NatID_filename,'NatIdGrid'); % read the NatID grid
    nc.NatIDVect = reshape(nc.NatIDGrid,[1 numel(nc.NatIDGrid)]); % as 1-D vect
    nc.land_point= ~isnan(nc.NatIDVect); % find land points in grid
else
    fprintf('WARNING (severe), file not found: %s\n',params.NatID_filename)
    fprintf('> might lead to (serious trouble)\n')
end

% deal with the (global) centroids
% --------------------------------

[~,fN]=fileparts(params.NatID_filename);
full_centroids_file=[centroids_file fN '.mat'];
params.centroids_file=full_centroids_file; % output only
if exist(full_centroids_file,'file')
    if verbose,fprintf('loading global isimip centroids from %s (delete this file to re-generate)\n',full_centroids_file);end
    load(full_centroids_file) % loads the struct centroids
else
    if isfield(nc,'land_point')
        fprintf('creating global isimip centroids file:\n')
        centroids.NatID_RegID=NatID_RegID;
        centroids.NatID=nc.NatIDVect(nc.land_point); % constrain to land and convert to 1D
        centroids.lon=vectlon(nc.land_point); % constrain to land and convert to 1D
        centroids.lat=vectlat(nc.land_point); % constrain to land and convert to 1D
        centroids.centroid_ID=1:length(centroids.lat); % define the GLOBAL centroid_ID
        
        n_centroids=length(centroids.lon);
        if params.distance_to_coast
            fprintf('> calculating distance to coast of %i land points, takes time...\n',n_centroids)
            centroids.distance2coast_km=climada_distance2coast_km(centroids.lon,centroids.lat);
        end
        
        fprintf('> adding coarse (1x1deg) ocean grid\n')
        for lon_i=-180:180 % add regular coarse grid
            for lat_i=-60:60
                centroids.lon(end+1)=lon_i+0.05; % to keep this points unique
                centroids.lat(end+1)=lat_i+0.05;
            end
        end
        centroids.centroid_ID=1:length(centroids.lon);
        centroids.centroid_ID(n_centroids+1:end)=centroids.centroid_ID(n_centroids+1:end)+(3e6-n_centroids); % to set the ocean-point centroids apart
        if params.distance_to_coast,centroids.distance2coast_km(n_centroids+1:length(centroids.lon))=0;end % ocean points
        
        fprintf('> saving global isimip centroids in %s\n',full_centroids_file);
        save(full_centroids_file,'centroids')
        
        if length(centroids.lon)>4000000 && params.distance_to_coast
            n_centroids=length(centroids.lon);centroids_full=centroids;
            pos=find(centroids.lat<60.06 & centroids.distance2coast_km<500);
            [centroids,untreated_fields]=climada_subarray(centroids,pos);
            if sum(strncmp('NatID',untreated_fields,5))>0,...
                    centroids.NatID = centroids.NatID(pos(pos<=length(centroids.NatID)));end
            [fP,fN,fE]=fileparts(full_centroids_file);
            full_centroids_file_red=[fP filesep fN '_red' fE];
            centroids.comment=sprintf('as %s, but latitude -60..60 and dense points only closer than 500km to coast, coarse grid also inland',fN);
            fprintf('> saving reduced (%i insted of %i centroids) global isimip centroids in %s\n',...
                n_centroids,length(pos),full_centroids_file_red);
            save(full_centroids_file_red,'centroids')
            centroids=centroids_full; centroids_full=[]; % keep going with full set
        end
        
    else
        fprintf('WARNING: unable to create proper global centroids (no NatID_file), using RAW method\n');
        centroids.lon=vectlon;
        centroids.lat=vectlat;
        centroids.centroid_ID=1:length(centroids.lat); % define the RAW GLOBAL centroid_ID
        full_centroids_file=[centroids_file fN 'RAW.mat'];
        fprintf('> saving RAW global isimip centroids in %s\n',full_centroids_file);
        save(full_centroids_file,'centroids');
    end
end % exist(centroids_file,'file')

% get template entity
entity=climada_entity_load(entity_template);

if strcmpi(ISO3,'ALL_IN_ONE')
    
    entity.assets.filename=[climada_global.entities_dir filesep 'GLB_' grid_resolution_str '_entity'];
    
    entity.assets.NatID_RegID=NatID_RegID;
    entity.assets.NatID=nc.NatIDVect(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lon=vectlon(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lat=vectlat(nc.land_point); % constrain to land and convert to 1D
    
    entity.assets.centroid_index=centroids.centroid_ID(1:length(entity.assets.lat));
    
    n_centroids=sum(sum(nc.land_point));
    if verbose,fprintf('> extracting %i land points at %i times (%i..%i) ...',n_centroids,n_times,time_val_yyyy(1),time_val_yyyy(end));end
    entity.assets.Values=zeros(n_times,n_centroids);
    if ~isempty(time_pop_yyyy),entity.assets.Population=zeros(n_times,n_centroids);end;pop_missing_count=0;
    if ~isempty(time_pop2_yyyy),entity.assets.Population2=zeros(length(time_pop2_yyyy),n_centroids);end;pop2_missing_count=0;
    for time_i=1:n_times
        temp_data=ncread(params.val_filename,val_variable_name,[1 1 time_i],[Inf Inf 1]); % only one time slab
        temp_data=temp_data.*nc.con_factor; % aply conversion factor
        temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
        entity.assets.Values(time_i,:)=temp_data(nc.land_point); % constrain to land
        
        if ~isempty(time_pop_yyyy)
            pop_time_i=find(time_pop_yyyy==time_val_yyyy(time_i));
            if length(pop_time_i)==1
                temp_data=ncread(params.pop_filename,pop_variable_name,[1 1 pop_time_i],[Inf Inf 1]); % only one time slab
                temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
                entity.assets.Population(time_i,:)=temp_data(nc.land_point); % constrain to land
            else
                pop_missing_count=pop_missing_count+1;
            end
        end % population
        
        if ~isempty(time_pop2_yyyy)
            pop2_time_i=find(time_pop2_yyyy==time_val_yyyy(time_i));
            if length(pop2_time_i)==1
                temp_data=ncread(params.pop2_filename,pop_variable_name,[1 1 pop2_time_i],[Inf Inf 1]); % only one time slab
                temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
                entity.assets.Population2(time_i,:)=temp_data(nc.land_point); % constrain to land
            else
                pop2_missing_count=pop2_missing_count+1;
            end
        end % population2
        
    end  % time_i
    if verbose,fprintf(' done\n');end
    
    if verbose && pop_missing_count>0, fprintf('WARNING: no population for        %2.2i years\n',pop_missing_count );end
    if verbose && pop2_missing_count>0,fprintf('WARNING: no second population for %2.2i years\n',pop2_missing_count);end

    entity.assets.isimip_comment='isimip whole world entity';
    
else
    
    if strcmpi(ISO3,'all')
        % process all
        n_NatIDs=length(NatID_RegID.NatID);
        if verbose,fprintf('processing %i countries (producing single country entity files)\n',n_NatIDs);end
        for iso3_i=1:n_NatIDs
            ISO3=NatID_RegID.ISO3{iso3_i};
            local_params.val_filename=val_filename;
            local_params.con_filename=con_filename;
            local_params.NatID_filename=NatID_filename;
            local_params.check_plot=check_plot;
            entity=isimip_gdp_entity(ISO3,local_params);
        end
        if verbose,fprintf('only last entity returned, see climada_entity_load\n');end
        return
    end
    
    if isempty(ISO3) % prompt for
        [selection] = listdlg('PromptString','Select country:',...
            'ListString',NatID_RegID.ISO3,'SelectionMode','multiple');
        pause(0.1)
        if ~isempty(selection)
            ISO3=NatID_RegID.ISO3(selection);
        else
            return % cancel pressed
        end
    end
    
    [~,~,iso3_pos]=intersect(ISO3,NatID_RegID.ISO3);
    NatID_RegID=climada_subarray(NatID_RegID,iso3_pos);
    
    if ~isempty(iso3_pos)
        NatID=NatID_RegID.NatID;
        entity.assets.NatID_RegID=NatID_RegID; % store the list
    end % ~isempty(iso3_pos)
    
    if iscell(ISO3) % more than one country
        ISO3=sort(ISO3); % sort to have alphabetical order
        ISO3_char=cell2mat(ISO3);
    else
        ISO3_char=ISO3;
    end
    
    NatID_pos=[]; % init
    if exist(params.NatID_filename,'file')
        if sum(abs(size(temp_data(:,:,1))-size(nc.NatIDGrid))) == 0 % grid sizes the same
            NatID_pos=ismember(nc.NatIDVect,NatID); % return all positions
            entity.assets.NatID=nc.NatIDVect(NatID_pos); % store the NatID grid
            n_centroids=sum(NatID_pos); % logical
            
            if verbose,fprintf('%i NatID grid cells within country %s (conversion min/max: %2.1f/%2.1f)\n',...
                    n_centroids,ISO3_char,min(nc.con_factor(NatID_pos)),max(nc.con_factor(NatID_pos)));end
        else
            fprintf('WARNING: grid sizes do not match, 2nd approach (shape)\n');
        end
    else
        fprintf('WARNING: NatIDGrid file not found, 2nd approach (shape)\n');
    end
    
    if isempty(NatID_pos) && ~iscell(ISO3) % 2nd try, using shapes
        
        % read admin0 (country) shape file
        admin0_shapes=climada_shaperead(admin0_shape_file);
        
        % figure the country shape
        shape_i=[];
        for shape_ii=1:length(admin0_shapes)
            if strcmpi(admin0_shapes(shape_ii).ADM0_A3,ISO3),shape_i=shape_ii;end
        end % shape_ii
        
        if ~isempty(shape_i)
            NatID_pos=climada_inpolygon(vectlon,vectlat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,0);
            n_centroids=sum(NatID_pos);
            if n_centroids>0
                if verbose,fprintf('%i grid cells within country %s shape\n',n_centroids,ISO3);end
            else
                fprintf('ERROR: no gridpoints within shape, aborted\n');
                return
            end
        else
            fprintf('ERROR: %s not in shape file, aborted\n',ISO3);
        end
    end
    
    if ~isempty(NatID_pos)
        
        entity.assets.filename=[climada_global.entities_dir filesep strtrim(ISO3_char) '_' grid_resolution_str '_entity'];
        entity.assets.admin0_ISO3=ISO3;
        
        entity.assets.lon=vectlon(NatID_pos);
        entity.assets.lat=vectlat(NatID_pos);
        
        entity.assets.centroid_index=centroids.centroid_ID(NatID_pos(nc.land_point));
        
        if verbose,fprintf('> extracting %i centroids at %i times (%i..%i) ...',n_centroids,n_times,time_val_yyyy(1),time_val_yyyy(end));end
        entity.assets.Values=zeros(n_times,n_centroids);
        if ~isempty(time_pop_yyyy),entity.assets.Population=zeros(n_times,n_centroids);end;pop_missing_count=0;
        if ~isempty(time_pop2_yyyy),entity.assets.Population2=zeros(length(time_pop2_yyyy),n_centroids);end;pop2_missing_count=0;
        for time_i=1:n_times
            temp_data=ncread(params.val_filename,val_variable_name,[1 1 time_i],[Inf Inf 1]); % only one time slab
            temp_data=temp_data.*nc.con_factor; % aply conversion factor
            temp_data=reshape(temp_data,[1 numel(temp_data)]); % convert to 1D
            entity.assets.Values(time_i,:)=temp_data(NatID_pos); % store
            
            if ~isempty(time_pop_yyyy)
                pop_time_i=find(time_pop_yyyy==time_val_yyyy(time_i));
                if length(pop_time_i)==1
                    temp_data=ncread(params.pop_filename,pop_variable_name,[1 1 pop_time_i],[Inf Inf 1]); % only one time slab
                    temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
                    entity.assets.Population(time_i,:)=temp_data(NatID_pos); % store
                else
                    pop_missing_count=pop_missing_count+1;
                end
            end % population
            
            if ~isempty(time_pop2_yyyy)
                pop2_time_i=find(time_pop2_yyyy==time_val_yyyy(time_i));
                if length(pop2_time_i)==1
                    temp_data=ncread(params.pop2_filename,pop_variable_name,[1 1 pop2_time_i],[Inf Inf 1]); % only one time slab
                    temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
                    entity.assets.Population2(time_i,:)=temp_data(NatID_pos); % store
                else
                    pop2_missing_count=pop2_missing_count+1;
                end
            end % population2
        
        end  % time_i
        if verbose,fprintf(' done\n');end
        
        if verbose && pop_missing_count>0, fprintf('WARNING: no population for        %2.2i years\n',pop_missing_count );end
        if verbose && pop2_missing_count>0,fprintf('WARNING: no second population for %2.2i years\n',pop2_missing_count);end
        
        entity.assets.isimip_comment=sprintf('isimip entity, created %s',datestr(now));
        
    end % ~isempty(NatID_pos)
    
end % ALL_IN_ONE

if isfield(entity.assets,'isimip_comment') % indicates we have an ok entity
    
    entity.assets.Values_yyyy=time_val_yyyy;
    
    % convert currency units
    entity.assets.Values=entity.assets.Values*Values_factor;
    
    % set active assets to first entry
    entity.assets.reference_year=entity.assets.Values_yyyy(1);
    entity.assets.Value=entity.assets.Values(1,:);
    
    % complete entity
    entity.assets.Cover      =entity.assets.Value;
    entity.assets.Deductible =entity.assets.Value*0;
    entity.assets.DamageFunID=entity.assets.Deductible+1;
    entity.assets.Category_ID=entity.assets.DamageFunID;
    entity.assets.Region_ID  =entity.assets.DamageFunID;
    
    if isfield(entity.assets,'Value_unit'),entity.assets=rmfield(entity.assets,'Value_unit');end
    if isfield(entity.assets,'hazard'),entity.assets=rmfield(entity.assets,'hazard');end
    
    entity.assets = climada_assets_complete(entity.assets);
        
    % add the source files
    entity.assets.centroids_file=params.centroids_file; % store the centroid origin
    entity.assets.isimip_gdp_entity_params=params; % to really keep track
    
    if ~isempty(params.hazard_file)
        hazard=climada_hazard_load(params.hazard_file);
        if ~isempty(hazard)
            entity.hazard.lon             = hazard.lon(entity.assets.centroid_index);
            entity.hazard.lat             = hazard.lat(entity.assets.centroid_index);
            entity.hazard.centroid_ID     = hazard.centroid_ID(entity.assets.centroid_index);
            entity.hazard.intensity       = hazard.intensity(:,entity.assets.centroid_index);
            entity.hazard.fraction        = hazard.fraction(:, entity.assets.centroid_index);
            entity.hazard.frequency       = hazard.frequency;
            entity.hazard.yyyy            = hazard.yyyy;
            entity.hazard.name            = hazard.name;
            entity.hazard.comment         = hazard.comment;
            entity.hazard.peril_ID        = hazard.peril_ID;
            entity.hazard.units           = hazard.units;
            entity.hazard.filename        = hazard.filename;
            entity.hazard.date            = hazard.date;
            entity.hazard.annotation_str  = hazard.annotation_str;
            % follow non-key fields, but easier to keep full hazard structure
            entity.hazard.reference_year  = hazard.reference_year;
            entity.hazard.orig_years      = hazard.orig_years;
            entity.hazard.orig_event_count= hazard.orig_event_count;
            entity.hazard.event_count     = hazard.event_count;
            entity.hazard.event_ID        = hazard.event_ID;
            entity.hazard.orig_event_flag = hazard.orig_event_flag;
            
            entity.assets.global_centroid_index=entity.assets.centroid_index; % the index into the full hazard set (backup)
            entity.assets.centroid_index  =1:length(entity.hazard.lon);

            fprintf('hazard joined from %s\n',params.hazard_file);
            
            if params.hazard_match
                % reduce to the years we have assets and hazard
                [common_yyyy,iasset]=intersect(entity.assets.Values_yyyy,entity.hazard.yyyy);
                fprintf('both assets and hazard for years %i..%i\n',min(common_yyyy),max(common_yyyy));
                hazard_pos=ismember(entity.hazard.yyyy,common_yyyy);
                entity.hazard.yyyy            = entity.hazard.yyyy(hazard_pos); % restrict to years
                entity.hazard.intensity       = entity.hazard.intensity(hazard_pos,:);
                entity.hazard.fraction        = entity.hazard.fraction(hazard_pos,:);
                entity.hazard.frequency       = entity.hazard.frequency(hazard_pos);
                entity.hazard.name            = entity.hazard.name(hazard_pos);
                entity.hazard.event_ID        = entity.hazard.event_ID(hazard_pos);
                entity.hazard.orig_event_flag = entity.hazard.orig_event_flag(hazard_pos);
                entity.hazard.event_count     = length(entity.hazard.event_ID);
                entity.hazard.orig_event_count = sum(entity.hazard.orig_event_flag);
                entity.hazard.orig_years      = length(unique(entity.hazard.yyyy));
                entity.assets.Values          = entity.assets.Values(iasset,:);
                if isfield(entity.assets,'Population'),entity.assets.Population=entity.assets.Population(iasset,:);end
                if isfield(entity.assets,'Population2'),entity.assets.Population2=entity.assets.Population2(iasset,:);end
                entity.assets.Values_yyyy     = entity.assets.Values_yyyy(iasset);
                if sum(abs(unique(entity.hazard.yyyy)-unique(entity.assets.Values_yyyy)))>0
                    fprintf('WARNING: match between assets and hazard might be incomplete\n');
                end
            end % params.hazard_match
            
            % a unique ID, a number for faster comparison, hence convert 1950166N14262 to 1950166.14262
            n_events=length(entity.hazard.name);
            entity.hazard.ID_no=zeros(1,n_events);
            for event_i=1:n_events
                entity.hazard.ID_no(event_i)=str2double(entity.hazard.name{event_i}(1:7))+str2double(entity.hazard.name{event_i}(9:end))/100000;
            end % event_i
            
        else
            fprintf('WARNING: hazard not joined (%s not found)\n',params.hazard_file);
        end
    end % ~isempty(params.hazard_file)
    
    if verbose,fprintf('> saving entity as %s\n',entity.assets.filename);end
    save(entity.assets.filename,'entity','-v7.3'); % -v7.3 for size...
    
    if params.check_plot
        climada_entity_plot(entity);
        % plot the sum of assets per year
        val_yearly_sum=sum(entity.assets.Values,2);
        [~,fN]=fileparts(entity.assets.filename);
        title_str=strrep(fN,'_entity','');
        figure('Name','Values checksum','Color',[1 1 1]),bar(entity.assets.Values_yyyy,val_yearly_sum,'FaceColor',[0 .5 .5],'EdgeColor',[1 1 1]);
        axis tight;title(title_str);xlabel('years');ylabel('yearly sum of assets')
    end
else
    entity=[];
    return
end

end % isimip_gdp_entity