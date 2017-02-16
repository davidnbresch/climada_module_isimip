function [entity,params]=isimip_gdp_entity(ISO3,params)
% climada isimip entity population gdp
% MODULE:
%   isimip
% NAME:
%   isimip_gdp_entity
% PURPOSE:
%   create the entity based on isimip GDP or population data.
%
%   Reads the GDP or population grid (per time, usually once per year) from
%   data_filename, multiplies with the conversion_factor (per gridcell)
%   from conversion_filename and stores into climada
%   entity.assets.Values(i,j) for year i and centroid j. The code maps to
%   the GLOBAL centroids, such that one can generate one (global) hazard
%   set and all single/multi country entities run properly (see
%   entity.assets.centroid_index on otuput below).
%
%   Please be PATIENT the first time you run this code, as it generates the
%   global reference grid, where calculation of distance to coast for all
%   land points does take (substantial) time. See PARAMETERS to switch
%   distance_to_coast off. In the first run, it also saves a centroids
%   file, in case more than 4 mio centroids, it also saves a reduced
%   verison (see code for details).
%
%   single/multi country mode (e.g. ISO3='DEU', ISO3={'DEU','FRA'}):
%    Checks the country ID(s) (NatId) on NatId_filename and takes all gridcells
%    within the requested country(s). If the resolution of the NatId does not
%    match the population data or if there is no NatId_filename provided,
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
% INPUTS:
%   ISO3: the ISO3 country code, e.g. 'DEU' or 'FRA' or {'DEU','FRA'} to
%       combine two countries in one entity (the entity still knows which
%       centroids belong to which country, see entity.assets.NatId)
%       if ='all': process all countries (be careful, test a few first),
%           create one single entity for each country (takes time and memory,
%           hence see also nexst option).
%       if ='ALL_IN_ONE', one entity with all countries is created
%           This is a bit a special case, as the entity is mainly to be
%           used within the ismip context in this case.
%       > promted for (to select from a list, also multiple) if not given
%       if ='params', just return default parameters, in entity, i.e. the
%           first output, already.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    data_filename: filename of the .nc file with the GDP or population values
%       ='0360as': use 0.1 degree default file (default, set in PARAMETERS)
%       ='0150as': use 2.5 minutes default file (set in PARAMETERS)
%       If only a filename (without path) is passed, isimip_data_dir (set
%       in PARAMETERS) is prepended
%    conversion_filename: filename of the .nc file with the conversion
%       factor from either GDP or population to assets (per gridpoint)
%       if data_filename is one of the short names, conversion_filename is
%       set to the corresponding file, if empty on input. Default set in
%       PARAMETERS. If only a filename (without path) is passed,
%       isimip_data_dir (set in PARAMETERS) is prepended.
%    NatId_filename: filename of the .nc file with the national grid IDs (to
%       assign grid cells to countries). Default set in PARAMETERS
%       if data_filename is one of the short names, NatId_filename is
%       set to the corresponding file, if empty on input
%       To avoid using NatId and just use the grid cells within a country
%       shape, set NatId_filename='ignore', but if ISO3='ALL_IN_ONE', NatId
%       is needed, hence ='ignore' does not work.
%       If only a filename (without path) is passed, isimip_data_dir (set
%       in PARAMETERS) is prepended .
%    check_plot: whether show a check plot (=1, default), or not (=0)
% % %    time_t0: to set the first year (as variable 'time' does not contain it)
% % %       default=1980, if time(1)<1800|>2100, otherwise time_t0=0
%    distance_to_coast: whether we calculate distant to coast (in km) of
%       call centroids (default=1), which speeds up later climada calculations
%       substantially (as coastal hazards need not to be evaluated at
%       inner-continental points). Set =0 in special cases, as initial
%       calculation might easily take 1h (we're talking about millions of
%       centroids... but since climada_distance2coast_km listens to
%       climada_global.parfor, set climada_global.parfor=1 for speedup).
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
%   PLUS the fields
%       entity.assets.Values(n_times,n_centroids) with the
%           variable as on the gdp file for each timestep at each centroid.
%           These values are NOT scaled (different from entity.assets.Value,
%           which is scaled by GDP and income group).
%       entity.assets.Values_yyyy: the year the Values(i,:) are valid for
%       entity.assets.NatId: the isimip country number for each centroid
%       entity.assets.ISO3_list: the list linking isimip country numbers
%           in ISO3_list(:,2) with ISO names in ISO3_list(:,1)
%       Note that entity.assets.centroid_index is already pointing to the
%           correct centroids in the corresponding global centroids file
%           (as retunred in params.centroid_file), but the hazard is not
%           set yet (instead of re-encoding, just define
%           entity.assets.hazard yourself).
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161017, initial
% David N. Bresch, david.bresch@gmail.com, 20161120, ALL_IN_ONE added
% David N. Bresch, david.bresch@gmail.com, 20161121, list of ISO3 allowed
% David N. Bresch, david.bresch@gmail.com, 20161123, checked for 2.5min
% David N. Bresch, david.bresch@gmail.com, 20170205, renamed to isimip_gdp_entity, overhaul, conversion_factor, params
% David N. Bresch, david.bresch@gmail.com, 20170209, 0150as files updated
% David N. Bresch, david.bresch@gmail.com, 20170215, ISO3 sorted, if a list
%-

entity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('ISO3','var'),   ISO3   = '';end
if ~exist('params','var'), params = struct;end

% check for some parameter fields we need
if ~isfield(params,'data_filename'),      params.data_filename='';end
if ~isfield(params,'conversion_filename'),params.conversion_filename='';end
if ~isfield(params,'NatId_filename'),     params.NatId_filename='';end
if ~isfield(params,'check_plot'),         params.check_plot=[];end
% if ~isfield(params,'time_t0'),            params.time_t0=[];end
if ~isfield(params,'distance_to_coast'),  params.distance_to_coast=[];end
if ~isfield(params,'centroids_file'),     params.centroids_file='';end % output only

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
% define default filenames for population or gdp, conversion and country number (NatId)
%data_filename0360as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0360as_yearly_zip.nc'];data_variable_name='var1';
data_filename0360as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0360as_remapnn_yearly.nc'];data_variable_name='gdp_grid';
conversion_filename0360as=[isimip_data_dir filesep 'GDP2Asset_converter_0360as_adv_1.nc']; % conversion_factor
%NatId_filename0360as=[isimip_data_dir filesep 'Nat_id_grid_0360as.nc']; % no tolerance arund borders
NatId_filename0360as=[isimip_data_dir filesep 'Nat_id_grid_0360as_adv_1.nc'];
%
%data_filename0150as =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0150as_yearly_zip.nc4']; % data_variable_name='var1';
data_filename0150as =[isimip_data_dir filesep 'gdp_1980-2100_SSP2_0150as_remapnn_yearly.nc']; % data_variable_name='var1';
conversion_filename0150as=[isimip_data_dir filesep 'GDP2Asset_converter_0150as.nc']; % conversion_factor
NatId_filename0150as=[isimip_data_dir filesep 'Nat_id_grid_0150as.nc'];
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
% if isempty(params.time_t0),           params.time_t0=[];end
if isempty(params.distance_to_coast), params.distance_to_coast=1;end

if strcmpi(ISO3,'params'),entity=params;return;end % special case, return the full params structure

if isempty(params.data_filename),params.data_filename='0360as';end

if strcmpi(params.data_filename,'0360as')
    params.data_filename=data_filename0360as;
    if isempty(params.conversion_filename),params.conversion_filename=conversion_filename0360as;end
    if isempty(params.NatId_filename),params.NatId_filename=NatId_filename0360as;end
elseif strcmpi(params.data_filename,'0150as')
    params.data_filename=data_filename0150as;
    if isempty(params.conversion_filename),params.conversion_filename=conversion_filename0150as;end
    if isempty(params.NatId_filename),params.NatId_filename=NatId_filename0150as;end
end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.data_filename,'file'),params.data_filename=[isimip_data_dir filesep params.data_filename];end
if ~exist(params.conversion_filename,'file'),params.conversion_filename=[isimip_data_dir filesep params.conversion_filename];end
if ~exist(params.NatId_filename,'file'),params.NatId_filename=[isimip_data_dir filesep params.NatId_filename];end

% read the population or GDP data
% read only one slide (timestep) to save memory, read slab-by-slab below
% -------------------------------

% first, figure the variable name for the data
nc.info=ncinfo(params.data_filename);

% figure value units
Values_factor=1;
try
    if strcmpi(nc.info.Variables(4).Attributes(2).Value,'GDP PPP 2005 billion $ US')
        Values_factor=1e9;
    end
catch
    if verbose,fprintf('Warning: Determining Value unit from nc failed\n');end
end % try
if verbose,fprintf('Values factor f=%g (climada Values=f * Values from nc)\n',Values_factor);end

for i=1:length(nc.info.Variables)
    if strcmpi(nc.info.Variables(i).Name,'var1'),data_variable_name='var1';end
    if strcmpi(nc.info.Variables(i).Name,'gdp_grid'),data_variable_name='gdp_grid';end
end % i

if verbose,fprintf('reading lon, lat, time and %s from %s ...',data_variable_name,params.data_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon       = ncread(params.data_filename,'lon')';
nc.lat       = ncread(params.data_filename,'lat')';
nc.time      = ncread(params.data_filename,'time')';
n_times      = length(nc.time);
if TEST_mode,n_times=min(3,n_times);nc.time=nc.time(1:n_times);end % TEST mode, read only few times
nc.var1      = ncread(params.data_filename,data_variable_name,[1 1 1],[Inf Inf 1]); % only one time slab
if verbose,fprintf(' done\n');end

if verbose,fprintf('reading conversion_factor from %s\n',params.conversion_filename);end
nc.conversion_factor = ncread(params.conversion_filename,'conversion_factor');
if sum(abs(size(nc.var1)-size(nc.conversion_factor)))>0
    fprintf('Error: sizes of %s and conversion_factor do not match, aborted\n',data_variable_name);
    return
else
    fprintf('  conversion_factor, global min/max = %2.2g %2.2g\n',min(min(nc.conversion_factor)),max(max(nc.conversion_factor)));
end

ISO3_list=isimip_ISO3_list; % get the mapping ISO3 - isimip country code

% create the grid
if verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon,gridlat] = meshgrid(nc.lon,nc.lat);
gridlon=gridlon';gridlat=gridlat'; % still as grid
vectlon=reshape(gridlon,[1 numel(gridlon)]); % as 1-D vect
vectlat=reshape(gridlat,[1 numel(gridlat)]);
if verbose,fprintf(' done\n');end

if exist(params.NatId_filename,'file')
    nc.NatIdGrid = ncread(params.NatId_filename,'NatIdGrid'); % read the NatId grid
    nc.NatIdVect = reshape(nc.NatIdGrid,[1 numel(nc.NatIdGrid)]); % as 1-D vect
    nc.land_point= ~isnan(nc.NatIdVect); % find land points in grid
else
    fprintf('Warning (severe), file not found: %s\n',params.NatId_filename)
    fprintf('> might lead to (serious trouble)\n')
end

% deal with the (global) centroids
% --------------------------------

[~,fN]=fileparts(params.NatId_filename);
full_centroids_file=[centroids_file fN '.mat'];
params.centroids_file=full_centroids_file; % output only
if exist(full_centroids_file,'file')
    if verbose,fprintf('loading global isimip centroids from %s (delete this file to re-generate)\n',full_centroids_file);end
    load(full_centroids_file) % loads the struct centroids
else
    if isfield(nc,'land_point')
        fprintf('creating global isimip centroids file:\n')
        centroids.ISO3_list=ISO3_list;
        centroids.NatId=nc.NatIdVect(nc.land_point); % constrain to land and convert to 1D
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
            if sum(strncmp('NatId',untreated_fields,5))>0,...
                    centroids.NatId = centroids.NatId(pos(pos<=length(centroids.NatId)));end
            [fP,fN,fE]=fileparts(full_centroids_file);
            full_centroids_file_red=[fP filesep fN '_red' fE];
            centroids.comment=sprintf('as %s, but latitude -60..60 and dense points only closer than 500km to coast, coarse grid also inland',fN);
            fprintf('> saving reduced (%i insted of %i centroids) global isimip centroids in %s\n',...
                n_centroids,length(pos),full_centroids_file_red);
            save(full_centroids_file_red,'centroids')
            centroids=centroids_full; centroids_full=[]; % keep going with full set
        end
        
    else
        fprintf('> Warning: unable to create proper global centroids (no NatId_file), using RAW method\n');
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
    
    entity.assets.filename=[climada_global.entities_dir filesep 'GLB_isimip_entity'];
    
    entity.assets.ISO3_list=ISO3_list;
    entity.assets.NatId=nc.NatIdVect(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lon=vectlon(nc.land_point); % constrain to land and convert to 1D
    entity.assets.lat=vectlat(nc.land_point); % constrain to land and convert to 1D
    
    entity.assets.centroid_index=centroids.centroid_ID(1:length(entity.assets.lat));
    
    n_centroids=sum(sum(nc.land_point));
    if verbose,fprintf('> extracting %i land points at %i times ...',n_centroids,n_times);end
    entity.assets.Values=zeros(n_times,n_centroids);
    for time_i=1:n_times
        temp_data=ncread(params.data_filename,data_variable_name,[1 1 time_i],[Inf Inf 1]); % only one time slab
        temp_data=temp_data.*nc.conversion_factor; % aply conversion factor
        temp_data=reshape(temp_data,[1 numel(temp_data)]); % as 1-D vect
        temp_data=temp_data(nc.land_point); % constrain to land
        entity.assets.Values(time_i,:)=temp_data; % store
    end  % time_i
    if verbose,fprintf(' done\n');end
    
    entity.assets.isimip_comment='isimip whole world entity, entity.assets.Value(s) not scaled';
    
else
    
    if strcmpi(ISO3,'all')
        % process all
        if verbose,fprintf('processing %i countries (producing single country entity files)\n',length(ISO3_list));end
        for iso3_i=1:length(ISO3_list)
            ISO3=ISO3_list(iso3_i,1);
            local_params.data_filename=data_filename;
            local_params.conversion_filename=conversion_filename;
            local_params.NatId_filename=NatId_filename;
            local_params.check_plot=check_plot;
            entity=isimip_gdp_entity(ISO3,local_params);
        end
        if verbose,fprintf('only last entity returned, see climada_entity_load\n');end
        return
    end
    
    if isempty(ISO3)
        % prompt for
        [selection] = listdlg('PromptString','Select country:',...
            'ListString',ISO3_list(:,1),'SelectionMode','multiple');
        pause(0.1)
        if ~isempty(selection)
            ISO3=ISO3_list(selection,1)';
        else
            return
        end
    end
    
    [~,~,iso3_pos]=intersect(ISO3,ISO3_list(:,1));
    
    if ~isempty(iso3_pos)
        NatId=cell2mat(ISO3_list(iso3_pos,2));
        entity.assets.ISO3_list=ISO3_list(iso3_pos,:); % store the list
    end % ~isempty(iso3_pos)
    
    if iscell(ISO3) % more than one country
        ISO3=sort(ISO3); % sort to have alphabetical order
        ISO3_char=cell2mat(ISO3);
    else
        ISO3_char=ISO3;
    end
    
    NatId_pos=[]; % init
    if exist(params.NatId_filename,'file')
        if sum(abs(size(nc.var1(:,:,1))-size(nc.NatIdGrid))) == 0 % grid sizes the same
            NatId_pos=ismember(nc.NatIdVect,NatId); % return all positions
            entity.assets.NatId=nc.NatIdVect(NatId_pos); % store the NatId grid
            n_centroids=sum(NatId_pos); % logical
            
            if verbose,fprintf(' %i NatId grid cells within country %s (conversion min/max: %2.1f/%2.1f)\n',...
                    n_centroids,ISO3_char,min(nc.conversion_factor(NatId_pos)),max(nc.conversion_factor(NatId_pos)));end
        else
            fprintf('Warning: grid sizes do not match, 2nd approach (shape)\n');
        end
    else
        fprintf('Warning: NatIdGrid file not found, 2nd approach (shape)\n');
    end
    
    if isempty(NatId_pos) && ~iscell(ISO3) % 2nd try, using shapes
        
        % read admin0 (country) shape file
        admin0_shapes=climada_shaperead(admin0_shape_file);
        
        % figure the country shape
        shape_i=[];
        for shape_ii=1:length(admin0_shapes)
            if strcmpi(admin0_shapes(shape_ii).ADM0_A3,ISO3),shape_i=shape_ii;end
        end % shape_ii
        
        if ~isempty(shape_i)
            NatId_pos=climada_inpolygon(vectlon,vectlat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,0);
            n_centroids=sum(NatId_pos);
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
    
    if ~isempty(NatId_pos)
        
        entity.assets.filename=[climada_global.entities_dir filesep ISO3_char '_entity'];
        entity.assets.admin0_ISO3=ISO3;
        
        entity.assets.lon=vectlon(NatId_pos);
        entity.assets.lat=vectlat(NatId_pos);
        
        entity.assets.centroid_index=centroids.centroid_ID(NatId_pos(nc.land_point));
        
        if verbose,fprintf(' extracting %i centroids at %i times ...',n_centroids,n_times);end
        entity.assets.Values=zeros(n_times,n_centroids);
        for time_i=1:n_times
            temp_data=ncread(params.data_filename,data_variable_name,[1 1 time_i],[Inf Inf 1]); % only one time slab
            temp_data=temp_data.*nc.conversion_factor; % aply conversion factor
            temp_data=reshape(temp_data,[1 numel(temp_data)]); % convert to 1D
            entity.assets.Values(time_i,:)=temp_data(NatId_pos); % store
        end  % time_i
        if verbose,fprintf(' done\n');end
        
        entity.assets.isimip_comment='isimip entity, entity.assets.Value and entity.assets.Values not modified (yet)';
        
    end % ~isempty(NatId_pos)
    
end % ALL_IN_ONE

if isfield(entity.assets,'isimip_comment') % indicates we have an ok entity
    
    entity.assets.Values_yyyy=zeros(1,n_times); % init
    for i=1:n_times
        entity.assets.Values_yyyy(i)=str2double(datestr(nc.time(i)+datenum(1950,1,1)+1,'yyyy'));
    end
    if entity.assets.Values_yyyy(1)<1860 || entity.assets.Values_yyyy(1)>2100 % strange content of nc.time
        entity.assets.Values_yyyy=1980+(1:n_times);
        fprintf('Warning: Values_yyyy hard wired to %i..%i\n',entity.assets.Values_yyyy(1),entity.assets.Values_yyyy(end))
    end
    if verbose,fprintf('timestamp %i .. %i\n',entity.assets.Values_yyyy(1),entity.assets.Values_yyyy(end));end
    
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
    entity.assets.centroids_file=params.centroids_file; % store the centroid origin
    
    entity.assets = climada_assets_complete(entity.assets);
    
    % add the source files
    entity.assets.data_filename =params.data_filename;
    entity.assets.NatId_filename=params.NatId_filename;
    
    if verbose,fprintf('saving entity as %s\n',entity.assets.filename);end
    save(entity.assets.filename,'entity','-v7.3'); % -v7.3 for size...
    
    if params.check_plot
        climada_entity_plot(entity);
        % plot the sum of assets per year
        value_yearly_sum=sum(entity.assets.Values,2);
        [~,fN]=fileparts(entity.assets.filename);
        title_str=strrep(fN,'_entity','');
        figure('Name','Values checksum','Color',[1 1 1]),bar(entity.assets.Values_yyyy,value_yearly_sum,'FaceColor',[0 .5 .5],'EdgeColor',[1 1 1]);
        axis tight;title(title_str);xlabel('years');ylabel('yearly sum of assets')
    end
else
    entity=[];
    return
end

end % isimip_gdp_entity