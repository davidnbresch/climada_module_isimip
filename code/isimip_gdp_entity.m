function [entity,params]=isimip_gdp_entity(ISO3,params,first_year,last_year,add_population)
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
%   Consider running entity=isimip_admin1_layer(entity) to add admin1
%   information to centroids.
%
%   Please be PATIENT the first time you run this code, as it generates the
%   global reference grid, where calculation of distance to coast for all
%   land points does take (substantial) time. See PARAMETERS to switch
%   distance_to_coast off. In the first run, it also saves a centroids
%   file, in case more than 4 mio centroids, it also saves a reduced
%   verison (see code for details).
%   HINT: Consider running isimip_admin1_layer to add admin1 information to centroids.
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
%    or: pkg install -params.verbose -forge -auto netcdf
%
%   next call: isimip...
% CALLING SEQUENCE:
%   entity=isimip_gdp_entity(ISO3,params)
% EXAMPLE:
%   entity=isimip_gdp_entity('DEU') % single country entity
%   entity=isimip_gdp_entity({'DEU','FRA','ITA'}) % multi country entity
%   entity=isimip_gdp_entity('all','0150as',1900,2020) % single country entity for each country worldwide in 0150as resolution, 1900..2020 (do NOT try this first ;-)
%   entity=isimip_gdp_entity('ALL_IN_ONE') % one global entity
%   params=isimip_gdp_entity('params') % get default parameters
%   % to check population for year 2030:
%   params.plot_population=1;params.year=2030;climada_entity_plot(entity,2,params)
%
%   % store the hazard also into the entity for special isimip use
%   params.hazard_file='GLB_0360as_TC_hist';
%   entity=isimip_gdp_entity('BRB',params)
% INPUTS:
%   ISO3: the ISO3 country code, e.g. 'DEU' or 'FRA' or {'DEU','FRA'} to
%       combine two countries in one entity (the entity still knows which
%       centroids belong to which country, see entity.assets.NatID)
%       if ='all': process all countries (be careful, test a few first),
%           create one single entity for each country.
%           Pro: keeps single entities handy (for e.g. 0150as). Consider
%           climada_global.parfor for substantial speedup.
%           Con: takes time to run and makes global calculations a bit
%           cumbersome, see also next option 'ALL_IN_ONE', at least for 0360as.
%       if ='ALL_IN_ONE', one entity with all countries is created
%           This is a bit a special case, as the entity is mainly to be
%           used within the ismip context.
%       > promted for (to select from a list, also multiple) if not given
%       if ='params', just return default parameters, in entity, i.e. the
%           first output, already.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    SPECIAL: if params is NOT a structure, it can be used to just define
%       the resolution, i.e. '0150as' or '0360as' (same effect as setting
%       params.val_filename='0360as' and then pass params, just easier).
%    grid_resolution: ='0360as' or ='0150as', this defines the resolution and
%       sets all files not passed via params accordingly. If the user passes
%       files via params (val_filename,pop_filename,pop2_filename...) it is
%       the user's repsonsibility to check for consistency. Often easiest
%       to just pass resolutoon as 2nd argument, see SPECIAL just above.
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
%    currency_unit: default=1, can be eg 1e6 or 1e9, in which case all
%       Values are divided by this unit before storing to the file, to store
%       hence millions or billions. Adds the field entity.assets.currency_unit
%       Note that other units than one on input (e.g. on the netCDF files)
%       are first converted to units of one, then currency_unit gets
%       applied.
%    verbose: if=1 (defualt), rather verbose, =2 for TEST mode (only few
%       timesteps processed)
%    entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%       entity_prefix='Try1' will result in Try1_DEU_0150as.mat
%   first_year: the first year to store, if not passed or empty (=[]), the
%       first year in the val (gdp) dataset, see there.
%   last_year: the last year to store, if not passed or empty (=[]), the
%       last year in the val (gdp) dataset, see there.
%   add_population: if =1, add population, =0 not (default)
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
%       entity.hazard.ID_no(event_i): ibtracsID as number, such as 1971275N10176 -> 1971275010176 (N???>0, S???>1)
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
% David N. Bresch, david.bresch@gmail.com, 20170304, pop2_* added
% David N. Bresch, david.bresch@gmail.com, 20170320, hazard.ID_no as integer
% David N. Bresch, david.bresch@gmail.com, 20170705, climada_global.save_file_version instead of hard-wired HDF5
% David N. Bresch, david.bresch@gmail.com, 20170706, currency_unit, gdp_1860-2100_0360as_yearly.nc and time_val_yyyy starting 1860
% David N. Bresch, david.bresch@gmail.com, 20180201, locate files
% David N. Bresch, david.bresch@gmail.com, 20180305, verbose=0 for 'all' and ready for latest 0150as file
% David N. Bresch, david.bresch@gmail.com, 20180316, time_val_yyyy for gdp calculated based on info in netCDF
% David N. Bresch, david.bresch@gmail.com, 20180316, Population_yyyy and Population2_yyyy added
% David N. Bresch, david.bresch@gmail.com, 20180317, parfor removed (might lead to parallel netCDF reading)
%-

entity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('ISO3','var'),   ISO3   = '';end
if ~exist('params','var'), params = struct;end
if ~exist('first_year','var'), first_year = [];end
if ~exist('last_year','var'),  last_year  = [];end
if ~exist('add_population','var'), add_population = [];end

if ischar(params) % special case, where we pass only the resolution in params
    grid_resolution=params;clear params;params = struct;
    params.grid_resolution=grid_resolution;clear grid_resolution;
end

% check for some parameter fields we need
if ~isfield(params,'grid_resolution'),    params.grid_resolution='';end
if ~isfield(params,'val_filename'),       params.val_filename='';end
if ~isfield(params,'con_filename'),       params.con_filename='';end
if ~isfield(params,'pop_filename'),       params.pop_filename='';end
if ~isfield(params,'pop2_filename'),      params.pop2_filename='';end
if ~isfield(params,'NatID_filename'),     params.NatID_filename='';end
if ~isfield(params,'check_plot'),         params.check_plot=[];end
if ~isfield(params,'distance_to_coast'),  params.distance_to_coast=[];end
if ~isfield(params,'hazard_file'),        params.hazard_file='';end
if ~isfield(params,'hazard_match'),       params.hazard_match=[];end
if ~isfield(params,'currency_unit'),      params.currency_unit=[];end
if ~isfield(params,'verbose'),            params.verbose=[];end
if ~isfield(params,'entity_prefix'),      params.entity_prefix='';end

% PARAMETERS
%
% to TEST the code, i.e. reads only a few (3) timesteps
TEST_mode=0; % default=0
%
if isempty(params.grid_resolution),params.grid_resolution='0360as';end % default='0360as'
%
if isempty(params.currency_unit)
    fprintf('WARNING: currency_unit set to 1.e9\n');
    params.currency_unit=1e9;
end
%
if isempty(params.verbose),params.verbose=1;end % default=1, to suppress output to stdout later
if params.verbose==2,params.verbose=1;TEST_mode=1;end % default=1, to suppress output to stdout later
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
%val_filename0360as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0360as_remapnn_yearly.nc']; % until 20170705
val_filename0360as =[isimip_data_dir filesep 'gdp_1860-2100_0360as_yearly.nc'];
%
val_variable_name='gdp_grid';
con_filename0360as=  [isimip_data_dir filesep 'GDP2Asset_converter_0360as_adv_1.nc'];
con_variable_name='conversion_factor';
pop_filename0360as = [isimip_data_dir filesep 'hyde_ssp2_1860-2100_0360as_yearly_zip.nc'];
pop2_filename0360as =[isimip_data_dir filesep 'hyde_ssp2_1860-2015_0360as_yearly_zip.nc4']; % pop_variable_name='gdp_grid';
pop_variable_name  = 'var1';
pop2_variable_name = 'var1';
NatID_filename0360as=[isimip_data_dir filesep 'NatID_grid_0360as_adv_1.nc'];
%gdp_first_year0360as  = 1860;
pop_first_year  = 1860;
pop2_first_year = 1860;
%
if isempty(add_population),add_population=0;end
%
%val_filename0150as =[isimip_data_dir filesep  'gdp_1860-2100_0150as_yearly.nc']; % val_variable_name='var1';
val_filename0150as =[isimip_data_dir filesep  'gdp_1850-2100_0150as_yearly.nc']; % val_variable_name='var1';
if ~exist(val_filename0150as,'file')
    if strfind(params.val_filename,'0150as'),fprintf('ERROR: wait for %s to be provided by Tobias\n',val_filename0150as);end
    val_filename0150as =[isimip_data_dir filesep  'gdp_1980-2100_SSP2_0150as_remapnn_yearly.nc'];                        % patch 20180201
    if strfind(params.val_filename,'0150as'),fprintf('--> PATCH, using %s for the time being\n',val_filename0150as);end % patch 20180201
    %gdp_first_year0150as  = 1950;
else
    %gdp_first_year0150as  = 1850; % 1850+(1:length(nc.time))-1
end
%val_filename0150as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0150as_remapnn_yearly.nc']; % val_variable_name='var1';

con_filename0150as=  [isimip_data_dir filesep 'GDP2Asset_converter_0150as.nc']; % con_factor
pop_filename0150as = [isimip_data_dir filesep 'hyde_ssp2_1860-2100_0150as_yearly_zip.nc4']; % pop_variable_name='gdp_grid';
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
% populate default parameters in params
if isempty(params.check_plot),        params.check_plot=0;end
if isempty(params.distance_to_coast), params.distance_to_coast=1;end
if isempty(params.hazard_match),      params.hazard_match=1;end
if isempty(params.currency_unit),     params.currency_unit=1;end

if strcmpi(params.grid_resolution,'0360as')
    %gdp_first_year=gdp_first_year0360as;
    if isempty(params.val_filename),  params.val_filename   = val_filename0360as;end
    if isempty(params.pop_filename),  params.pop_filename   = pop_filename0360as;end
    if isempty(params.pop2_filename), params.pop2_filename  = pop2_filename0360as;end
    if isempty(params.con_filename),  params.con_filename   = con_filename0360as;end
    if isempty(params.NatID_filename),params.NatID_filename = NatID_filename0360as;end
elseif strcmpi(params.grid_resolution,'0150as')
    %gdp_first_year=gdp_first_year0150as;
    if isempty(params.val_filename),  params.val_filename   = val_filename0150as;end
    if isempty(params.pop_filename),  params.pop_filename   = pop_filename0150as;end
    if isempty(params.pop2_filename), params.pop2_filename  = pop2_filename0150as;end
    if isempty(params.con_filename),  params.con_filename   = con_filename0150as;end
    if isempty(params.NatID_filename),params.NatID_filename = NatID_filename0150as;end
end

if strcmpi(ISO3,'params'),entity=params;return;end % special case, return the full params structure

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.val_filename,'file'),  params.val_filename   =[isimip_data_dir filesep params.val_filename];end
if ~exist(params.pop_filename,'file'),  params.pop_filename   =[isimip_data_dir filesep params.pop_filename];end
if ~exist(params.pop2_filename,'file'), params.pop2_filename  =[isimip_data_dir filesep params.pop2_filename];end
if ~exist(params.con_filename,'file'),  params.con_filename   =[isimip_data_dir filesep params.con_filename];end
if ~exist(params.NatID_filename,'file'),params.NatID_filename =[isimip_data_dir filesep params.NatID_filename];end

if ~add_population
    params.pop_filename='';
    params.pop2_filename='';
end

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
    if params.verbose,fprintf('WARNING: Determining Value unit from nc failed\n');end
end % try
if params.verbose,fprintf('Values factor f=%g (climada asset value = f * original value from file)\n',Values_factor);end

for i=1:length(nc.info.Variables)
    if strcmpi(nc.info.Variables(i).Name,'var1'),val_variable_name='var1';end
    if strcmpi(nc.info.Variables(i).Name,'gdp_grid'),val_variable_name='gdp_grid';end
end % i

if params.verbose,fprintf('reading lon, lat, time and %s from %s ...',val_variable_name,params.val_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon       = ncread(params.val_filename,'lon')';
nc.lat       = ncread(params.val_filename,'lat')';
nc.time      = ncread(params.val_filename,'time')';
n_times      = length(nc.time);
nc.time_units = ncreadatt(params.val_filename,'time','units');
temp_data    = ncread(params.val_filename,val_variable_name,[1 1 1],[Inf Inf 1]); % only one time slab
if params.verbose,fprintf(' done\n');end

% Get info about NetCDF time axis (used later)
nc.time_units = strsplit(nc.time_units, ' ');
nc.time_orig  = cell2mat(strcat(nc.time_units(3), {' '}, nc.time_units(4)));
nc.time_units = cell2mat(nc.time_units(1));

% create time axis
if strcmp(nc.time_units, 'years')
    mm=ones(1,n_times);
    dd=ones(1,n_times);
    time_val_yyyy=str2double(nc.time_orig(1:4))+nc.time';
    time_datenum=datenum(time_val_yyyy,mm,dd);
else
    if strcmp(nc.time_units, 'seconds')
        time_ratio = 3600*24;
    elseif strcmp(nc.time_units, 'minutes')
        time_ratio = 60*24;
    elseif strcmp(nc.time_units, 'hours')
        time_ratio = 24;
    elseif strcmp(nc.time_units, 'days')
        time_ratio = 1;
    else
        warning('** the time units in the NetCDF file are unexpected - assuming one every 366 days (roughly one year) ending in the reference year **')
        time_ratio=1/366;
        nc.time=(1:n_times)-1;
    end
    time_datenum = double(nc.time/time_ratio+datenum(nc.time_orig));
end

time_val_yyyy=str2num(datestr(time_datenum, 'yyyy'));

if min(time_val_yyyy)<1850 || max(time_val_yyyy)>2100 % strange content of nc.time
    fprintf('SEVERY WARNING: strange years %i..%i\n',min(time_val_yyyy),max(time_val_yyyy))
end

% restrict time
if ~isempty(first_year)
    pos=find(time_val_yyyy>=first_year);
    time_val_yyyy=time_val_yyyy(pos);
    time_val_start=pos(1);
else
    time_val_start=1;
end

if ~isempty(last_year)
    pos=find(time_val_yyyy<=last_year);
    time_val_yyyy=time_val_yyyy(pos);
    time_val_end=pos(end);
    time_val_end=time_val_start+time_val_end-1; % convert to absolute
else
    time_val_end=length(time_val_yyyy);
end

if TEST_mode % TEST mode, read only few times
    time_val_end=time_val_start+3;
    time_val_yyyy=time_val_yyyy(time_val_start:time_val_end);
end
n_times=time_val_end-time_val_start+1;

if params.verbose,fprintf('reading %s from %s',con_variable_name,params.con_filename);end
nc.con_factor = ncread(params.con_filename,con_variable_name);
if sum(abs(size(temp_data)-size(nc.con_factor)))>0
    if params.verbose,fprintf('\nError: sizes of %s and %s do not match, aborted\n',val_variable_name,con_variable_name);end
    return
else
    if params.verbose,fprintf(', min/max: %2.2g/%2.2g\n',min(min(nc.con_factor)),max(max(nc.con_factor)));end
end

if exist(params.pop_filename,'file')
    if params.verbose,fprintf('reading population from %s\n',params.pop_filename);end
    nc.time_pop=ncread(params.pop_filename,'time')';
    
    % currently hard-wired (as not properly defined in netCDF (there time just = 1..n)
    time_pop_yyyy=pop_first_year-1+(1:length(nc.time_pop));
    if params.verbose,fprintf('WARNING: population years hard wired to %i..%i\n',time_pop_yyyy(1),time_pop_yyyy(end));end
else
    time_pop_yyyy=[];
end

if exist(params.pop2_filename,'file')
    if params.verbose,fprintf('reading second population from %s\n',params.pop2_filename);end
    nc.time_pop2=ncread(params.pop2_filename,'time')';
    
    % currently hard-wired (as not properly defined in netCDF (there time just = 1..n)
    time_pop2_yyyy=pop2_first_year-1+(1:length(nc.time_pop2));
    if params.verbose,fprintf('WARNING: second population years hard wired to %i..%i\n',time_pop2_yyyy(1),time_pop2_yyyy(end));end
else
    time_pop2_yyyy=[];
end

NatID_RegID=isimip_NatID_RegID; % get the mapping ISO3 - isimip country code

% create the grid
if params.verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon,gridlat] = meshgrid(nc.lon,nc.lat);
gridlon=gridlon';gridlat=gridlat'; % still as grid
vectlon=reshape(gridlon,[1 numel(gridlon)]); % as 1-D vect
vectlat=reshape(gridlat,[1 numel(gridlat)]);
if params.verbose,fprintf(' done\n');end

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
    if params.verbose,fprintf('loading global isimip centroids from %s (delete this file to re-generate)\n',full_centroids_file);end
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
        save(full_centroids_file,'centroids',climada_global.save_file_version) % hdf5
        
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
            save(full_centroids_file_red,'centroids',climada_global.save_file_version) % hdf5
            centroids=centroids_full; clear centroids_full % keep going with full set
        end
        
    else
        fprintf('WARNING: unable to create proper global centroids (no NatID_file), using RAW method\n');
        centroids.lon=vectlon;
        centroids.lat=vectlat;
        centroids.centroid_ID=1:length(centroids.lat); % define the RAW GLOBAL centroid_ID
        full_centroids_file=[centroids_file fN 'RAW.mat'];
        fprintf('> saving RAW global isimip centroids in %s\n',full_centroids_file);
        save(full_centroids_file,'centroids',climada_global.save_file_version);
    end
end % exist(centroids_file,'file')

% get template entity
entity=climada_entity_load(entity_template);

% follow all options to request single, multiple or all countries

if isempty(ISO3) % prompt for country/-ies
    [selection] = listdlg('PromptString','Select country:',...
        'ListString',NatID_RegID.ISO3,'SelectionMode','multiple');
    pause(0.1)
    if ~isempty(selection)
        ISO3=NatID_RegID.ISO3(selection);
    else
        return % cancel pressed
    end
end

if strcmpi(ISO3,'ALL_IN_ONE')
    NatID_pos=[]; % no single country, one global entity
elseif strcmpi(ISO3,'all')
    % process all, call recursively
    n_NatIDs=length(NatID_RegID.NatID);
    params.verbose           = 0;
    params.check_plot        = 0;
    if TEST_mode,n_NatIDs=5;end % for TESTs only
    fprintf('processing %i countries (producing single country entity files)\n\n',n_NatIDs);
    for iso3_i=1:n_NatIDs
        ISO3=NatID_RegID.ISO3{iso3_i};
        fprintf('- %s: (%s)\n',ISO3,datestr(now))
        entity=isimip_gdp_entity(ISO3,params,first_year,last_year,add_population);
    end % iso3_i
    if params.verbose,fprintf('only last entity returned, see climada_entity_load\n');end
    return
else
    [~,~,iso3_pos]=intersect(ISO3,NatID_RegID.ISO3);
    NatID_RegID=climada_subarray(NatID_RegID,iso3_pos);
    
    if ~isempty(iso3_pos)
        NatID=NatID_RegID.NatID;
        entity.assets.NatID_RegID=NatID_RegID; % store the list
    else
        fprintf('ERROR %s not recognized\n',char(ISO3));
        return
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
            
            if params.verbose,fprintf('%i NatID grid cells within country %s (conversion min/max: %2.1f/%2.1f)\n',...
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
                if params.verbose,fprintf('%i grid cells within country %s shape\n',n_centroids,char(ISO3));end
            else
                fprintf('ERROR: no gridpoints within shape, aborted\n');
                return
            end
        else
            fprintf('ERROR: %s not in shape file, aborted\n',char(ISO3));
            return
        end
    end
end % to check all options for ISO3

if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end

if isempty(NatID_pos) % no single country, one global entity
    entity.assets.filename=[climada_global.entities_dir filesep params.entity_prefix 'GLB_' params.grid_resolution '_entity'];
    entity.assets.NatID_RegID=NatID_RegID;
    entity.assets.NatID=nc.NatIDVect(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lon=vectlon(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lat=vectlat(nc.land_point); % constrain to land and convert to 1D
    entity.assets.centroid_index=centroids.centroid_ID(1:length(entity.assets.lat));
    n_centroids=sum(sum(nc.land_point));
    if params.verbose,fprintf('> extracting %i land points at %i times (%i..%i) ... ',n_centroids,n_times,time_val_yyyy(1),time_val_yyyy(end));end
else % single country
    entity.assets.filename=[climada_global.entities_dir filesep params.entity_prefix strtrim(ISO3_char) '_' params.grid_resolution '_entity'];
    entity.assets.admin0_ISO3=ISO3;
    entity.assets.lon=vectlon(NatID_pos);
    entity.assets.lat=vectlat(NatID_pos);
    entity.assets.centroid_index=centroids.centroid_ID(NatID_pos(nc.land_point));
    n_centroids=sum(NatID_pos);
    if params.verbose,fprintf('> extracting %i centroids at %i times (%i..%i) ... ',  n_centroids,n_times,time_val_yyyy(1),time_val_yyyy(end));end
end

entity.assets.Values=zeros(n_times,n_centroids);
if ~isempty(time_pop_yyyy), entity.assets.Population =zeros(n_times,n_centroids);end;pop_missing_count=0;
if ~isempty(time_pop2_yyyy),entity.assets.Population2=zeros(1,      n_centroids);end;pop2_missing_count=0;

% find rectangle around country (doe snot work across dateline)
dlon=max(abs(diff(sort(nc.lon))));
dlat=max(abs(diff(sort(nc.lat))));
lonmin=min(entity.assets.lon)-dlon/2;
lonmax=max(entity.assets.lon)+dlon/2;
latmin=min(entity.assets.lat)-dlat/2;
latmax=max(entity.assets.lat)+dlat/2;
lonpos=find(nc.lon>=lonmin & nc.lon<=lonmax);
latpos=find(nc.lat>=latmin & nc.lat<=latmax);

temp_data_zeros=nc.NatIDGrid*0; % init an empty full slab

entity.assets.Values_yyyy=time_val_yyyy';

for time_i=time_val_start:time_val_end
    time_i_index=time_i-time_val_start+1; % index in output
    temp_data=temp_data_zeros; % init
    temp_data(lonpos,latpos)=ncread(params.val_filename,val_variable_name,...
        [lonpos(1) latpos(1) time_i],[length(lonpos) length(latpos) 1]); % only one time sub-slab
    temp_data=temp_data.*nc.con_factor; % aply conversion factor
    temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
    if isempty(NatID_pos)
        entity.assets.Values(time_i_index,:) = temp_data(nc.land_point); % constrain to land
    else
        entity.assets.Values(time_i_index,:) = temp_data(NatID_pos); % store
    end
    
    if ~isempty(time_pop_yyyy)
        pop_time_i=find(time_pop_yyyy==time_val_yyyy(time_i_index));
        if length(pop_time_i)==1
            entity.assets.Population_yyyy(time_i_index)=time_val_yyyy(time_i_index);
            temp_data=temp_data_zeros; % init
            temp_data(lonpos,latpos)=ncread(params.pop_filename,pop_variable_name,...
                [lonpos(1) latpos(1) pop_time_i],[length(lonpos) length(latpos) 1]); % only one time sub-slab
            temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
            if isempty(NatID_pos)
                entity.assets.Population(time_i_index,:)=temp_data(nc.land_point); % constrain to land
            else
                entity.assets.Population(time_i_index,:)=temp_data(NatID_pos); % store
            end
        else
            pop_missing_count=pop_missing_count+1;
        end
    end % population
    
    if ~isempty(time_pop2_yyyy)
        pop2_time_i=find(time_pop2_yyyy==time_val_yyyy(time_i_index));
        if length(pop2_time_i)==1
            entity.assets.Population2_yyyy(time_i_index)=time_val_yyyy(time_i_index);
            temp_data=temp_data_zeros; % init
            temp_data(lonpos,latpos)=ncread(params.pop2_filename,pop2_variable_name,...
                [lonpos(1) latpos(1) pop2_time_i],[length(lonpos) length(latpos) 1]); % only one time slab
            temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
            if isempty(NatID_pos)
                entity.assets.Population2(time_i_index,:)=temp_data(nc.land_point); % constrain to land
            else
                entity.assets.Population2(time_i_index,:)=temp_data(NatID_pos); % constrain to land
            end
        else
            pop2_missing_count=pop2_missing_count+1;
        end
    end % population2
    
end  % time_i
if params.verbose,fprintf('done\n');end

if params.verbose && pop_missing_count>0, fprintf('WARNING: no population for        %2.2i years\n',pop_missing_count );end
if params.verbose && pop2_missing_count>0,fprintf('WARNING: no second population for %2.2i years\n',pop2_missing_count);end

entity.assets.isimip_comment=sprintf('isimip entity, created %s',datestr(now));

% just an INFO
if params.verbose,fprintf('HINT: consider running entity=isimip_admin1_layer(entity) to add admin1 information to centroids\n');end

if isfield(entity.assets,'isimip_comment') % indicates we have an ok entity
    
    % convert currency units to ones
    entity.assets.Values=entity.assets.Values*Values_factor;
    
    % convert currency units to user requested unit
    if abs(params.currency_unit-1.0)>0
        entity.assets.Values=entity.assets.Values/params.currency_unit;
        entity.assets.currency_unit=params.currency_unit;
        if params.verbose,fprintf('NOTE: Value unit set to %g\n',params.currency_unit);end
    else
        entity.assets.currency_unit=1.0; % default in climada all in units of 1
    end
    
    % set active assets to first entry (default)
    entity.assets.Value=entity.assets.Values(1,:);
    entity.assets.reference_year=entity.assets.Values_yyyy(1);
    
    % if year 2005 exists, set default Value to 2005, as this is the
    % reference year for inflation
    pos=find(entity.assets.Values_yyyy==2005);
    if ~isempty(pos)
        entity.assets.Value=entity.assets.Values(pos,:);
        entity.assets.reference_year=entity.assets.Values_yyyy(pos);
    end
    
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
            
            % a unique ID, a number for faster comparison, hence convert 1950166N14262 to 1950166014262
            n_events=length(entity.hazard.name);
            entity.hazard.ID_no=zeros(1,n_events);
            for event_i=1:n_events
                %entity.hazard.ID_no(event_i)=str2double(entity.hazard.name{event_i}(1:7))+str2double(entity.hazard.name{event_i}(9:end))/100000; % 1950166N14262 to 1950166014262 (N???>0, S???>1)
                entity.hazard.ID_no(event_i)=str2double(strrep(strrep(entity.hazard.name{event_i},'S','1'),'N','0')); % N->0, S->1
            end % event_i
            
        else
            fprintf('WARNING: hazard not joined (%s not found)\n',params.hazard_file);
        end
    end % ~isempty(params.hazard_file)
    
    % fix the rounding errors in the centroid coordinates
    if strfind(params.val_filename,       '0360as')
        entity = fix_coords_isimip(entity,'0360as');
    elseif strfind(params.val_filename,   '0150as')
        entity = fix_coords_isimip(entity,'0150as');
    end
    
    fprintf('> saving entity as %s\n',entity.assets.filename);
    save(entity.assets.filename,'entity',climada_global.save_file_version); % HDF5
    
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