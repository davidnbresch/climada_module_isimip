function entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot,time_t0)
% climada isimip entity population
% MODULE:
%   isimip
% NAME:
%   isimip_ssp2_entity
% PURPOSE:
%   create the entity based on isimip ssp (population) data
%
%   single country mode (e.g. ISO3='DEU'):
%   Checks the country ID (NatId) on NatId_filename and takes all gridcells
%   within the requested country. If the resolution of the NatId does not
%   match the population data or if there is no NatId_filename provided,
%   the code uses the climada country shape files to select the gridcells
%   within the country (the code notifies to stdout).
%   The asset values are then scaled by GDP and the income group, see
%   climada_entity_value_GDP_adjust_one.
%
%   all country mode (ISO3='ALL_IN_ONE'):
%   same as for single country, but create one global entity with
%   additional fields as described in OUTPUTS below.
%
%   Octave: please install the netCDF package first:
%   pkg install -forge netcdf -local -auto
%   or: pkg install -verbose -forge -auto netcdf
%
%   next call: isimip...
% CALLING SEQUENCE:
%   entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot,time_t0)
% EXAMPLE:
%   entity=isimip_ssp2_entity('DEU') % single country entity
%   entity=isimip_ssp2_entity('ALL_IN_ONE') % one global entity
% INPUTS:
%   ISO3: the ISO3 country code, e.g. 'DEU' or 'FRA' or {'DEU','FRA'} to
%       combine two countries in one entity (the entity still knows which
%       centroids belong to which country, see entity.assets.NatId)
%       if ='all': process all countries (be careful, test a few first)
%       if ='ALL_IN_ONE', one entity with all countries is created
%           This is a bit a special case, as the entity is mainly to be
%           used within the ismip context in this case.
%       > promted for (to select from a list, also multiple) if not given
% OPTIONAL INPUT PARAMETERS:
%   ssp2_filename: filename of the .nc file with the population values
%       ='01deg': use 0.1 degree default file (set in PARAMETERS)
%       ='25min': use 2.5 minutes default file (set in PARAMETERS)
%   NatId_filename: filename of the .nc file with the national grid IDs (to
%       assign grid cells to countries). Default set in PARAMETERS
%       ='01deg': use 0.1 degree default file (set in PARAMETERS)
%       ='25min': use 2.5 minutes default file (set in PARAMETERS)
%       To avoid using NatId and just use the grid cells within a country
%       shape, set NatId_filename='ignore'.
%   check_plot: whether show a check plot (=1, default), or not (=0)
%   time_t0: to set the first year (as variable 'time' does not contain it)
%       default=1860, if time(1)<1800, otherwise time_t0=0
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
%   PLUS the fields
%       entity.assets.Values(n_times,n_centroids) with the
%           variable as on the ssp2 file for each timestep at each centroid.
%           These values are NOT scaled (different from entity.assets.Value,
%           which is scaled by GDP and income group).
%       entity.assets.Values_yyyy: the year the Values(i,:) are valid for
%       entity.assets.NatId: the isimip country number for each centroid
%       entity.assets.ISO3_list: the list linking isimip country numbers
%           in ISO3_list(:,2) with ISO names in ISO3_list(:,1)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161017, initial
% David N. Bresch, david.bresch@gmail.com, 20161120, ALL_IN_ONE added
% David N. Bresch, david.bresch@gmail.com, 20161121, list of ISO3 allowed
%-

entity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('ISO3','var'),          ISO3           = '';end
if ~exist('ssp2_filename','var'), ssp2_filename  = '';end
if ~exist('NatId_filename','var'),NatId_filename = '';end
if ~exist('check_plot','var'),    check_plot     =  0;end
if ~exist('time_t0','var'),       time_t0        =  [];end

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
% define default filenames for population (ssp2) and country number (NatId)
ssp2_filename01deg =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0.1deg_yearly_zip.nc'];
NatId_filename01deg=[isimip_data_dir filesep 'Nat_id_grid_0.1deg.nc'];
ssp2_filename25min =[isimip_data_dir filesep 'hyde_ssp2_1860-2100_2.5min_yearly_zip.nc'];
NatId_filename25min=[isimip_data_dir filesep 'Nat_id_grid_2.5min.nc'];
%
% the entity template to populate a defaiult entity
entity_template=[climada_global.entities_dir filesep 'entity_template'];
%
% admin0 shape file (for fallback selection option):
admin0_shape_file=climada_global.map_border_file;
%
% to TEST the code, i.e. reads only a few timesteps
TEST_mode=0; % default=0


if isempty(ssp2_filename),ssp2_filename='01deg';end
if isempty(NatId_filename),NatId_filename='01deg';end

if strcmpi(ssp2_filename,'01deg')
    ssp2_filename=ssp2_filename01deg;
    NatId_filename=NatId_filename01deg;
elseif strcmpi(ssp2_filename,'25min')
    ssp2_filename=ssp2_filename25min;
    NatId_filename=NatId_filename25min;
end

% read the population data
if verbose,fprintf('reading lon, lat, time and var1 from %s ...',ssp2_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon       = ncread(ssp2_filename,'lon');
nc.lat       = ncread(ssp2_filename,'lat');
nc.time      = ncread(ssp2_filename,'time');
n_times      = length(nc.time);
if TEST_mode,n_times=min(3,n_times);end % TEST mode, read only few times
nc.var1      = ncread(ssp2_filename,'var1',[1 1 1],[Inf Inf 1]); % only one time slab
if verbose,fprintf(' done\n');end

% create the grid
if verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon0,gridlat0] = meshgrid(nc.lon,nc.lat);
gridlon0=gridlon0';
gridlat0=gridlat0';

% get template entity
entity=climada_entity_load(entity_template);

ISO3_list=isimip_ISO3_list; % get the mapping ISO3 - isimip country code

if strcmpi(ISO3,'ALL_IN_ONE')
    
    % read the NatId grid
    nc.NatIdGrid = ncread(NatId_filename,'NatIdGrid');
    
    entity.assets.filename=[climada_global.entities_dir filesep 'isimip_entity'];
    
    land_point=~isnan(nc.NatIdGrid); % find land points
    entity.assets.NatId=nc.NatIdGrid(land_point)'; % constrain to land and convert to 1D
    entity.assets.ISO3_list=ISO3_list;
    entity.assets.lon=gridlon0(land_point);
    entity.assets.lat=gridlat0(land_point);
    
    n_centroids=sum(sum(land_point));
    if verbose,fprintf(' extracting %i land points at %i times ...',n_centroids,n_times);end
    entity.assets.Values=zeros(n_times,n_centroids);
    for time_i=1:n_times
        temp_var1=ncread(ssp2_filename,'var1',[1 1 time_i],[Inf Inf 1]); % only one time slab
        temp_var1=temp_var1(land_point); % constrain to land and convert to 1D
        entity.assets.Values(time_i,:)=temp_var1;
    end  % time_i
    if verbose,fprintf(' done\n');end
    
    entity.assets.isimip_comment='isimip whole world entity, entity.assets.Value(s) not scaled';
    
else
    
    if strcmpi(ISO3,'all')
        % process all
        if verbose,fprintf('processing %i entities\n',length(ISO3_list));end
        for iso3_i=1:length(ISO3_list)
            ISO3=ISO3_list(iso3_i,1);
            entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot);
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
    
    gridlon=reshape(gridlon0,[1 numel(gridlon0)]);
    gridlat=reshape(gridlat0,[1 numel(gridlat0)]);
    %if verbose,fprintf(' done\n');end
        
    NatId_pos=[]; % init
    if exist(NatId_filename,'file')
        % read the NatId grid
        nc.NatIdGrid = ncread(NatId_filename,'NatIdGrid');
        
        if sum(abs(size(nc.var1(:,:,1))-size(nc.NatIdGrid))) == 0 % grid sizes the same
            NatIdGrid=reshape(nc.NatIdGrid,[1 numel(nc.NatIdGrid)]);
            NatId_pos=ismember(NatIdGrid,NatId); % return all positions
            entity.assets.NatId=NatIdGrid(NatId_pos); % store the NatId grid
            n_centroids=sum(NatId_pos); % logical

            if verbose,fprintf(' %i NatId grid cells within country %s\n',n_centroids,cell2mat(ISO3));end
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
            NatId_pos=climada_inpolygon(gridlon,gridlat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,0);
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
        
        entity.assets.filename=[climada_global.entities_dir filesep cell2mat(ISO3) '_entity'];
        entity.assets.admin0_ISO3=ISO3;
        
        entity.assets.lon=gridlon(NatId_pos);
        entity.assets.lat=gridlat(NatId_pos);
        
        if verbose,fprintf(' extracting %i centroids at %i times ...',n_centroids,n_times);end
        entity.assets.Values=zeros(n_times,n_centroids);
        for time_i=1:n_times
            temp_var1=ncread(ssp2_filename,'var1',[1 1 time_i],[Inf Inf 1]); % only one time slab
            temp_var1=reshape(temp_var1,[1 numel(temp_var1)]); % convert to 1D
            entity.assets.Values(time_i,:)=temp_var1(NatId_pos);
        end  % time_i
        if verbose,fprintf(' done\n');end
        
        % scale up asset values based on a country's estimated total asset value
        %%entity=climada_entity_value_GDP_adjust_one(entity,verbose); %
        %%needs more sophisticated treatment, since GDP for each year etc.
        entity.assets.isimip_comment='isimip entity, entity.assets.Value and entity.assets.Values not modified (yet)';
        
    end % ~isempty(NatId_pos)
    
end % ALL_IN_ONE

if isfield(entity.assets,'isimip_comment') % indicates we have an ok entity
    
    if isempty(time_t0) && nc.time(1)<1800
        time_t0=1860;
    else
        time_t0=0;
    end
    entity.assets.Values_yyyy=nc.time(1:n_times)+time_t0;
    entity.assets.Values_yyyy=entity.assets.Values_yyyy'; % to have 1 x n_times
    if verbose,fprintf('timestamp %i .. %i\n',entity.assets.Values_yyyy(1),entity.assets.Values_yyyy(end));end
    
    % set active assets to first entry
    entity.assets.reference_year=entity.assets.Values_yyyy(1);
    entity.assets.Value=entity.assets.Values(1,:);
    
    % complete entity
    entity.assets.Cover      =entity.assets.Value;
    entity.assets.Deductible =entity.assets.Value*0;
    entity.assets.DamageFunID=entity.assets.Deductible+1;
    entity.assets.Category_ID=entity.assets.DamageFunID;
    entity.assets.Region_ID  =entity.assets.DamageFunID;
    
    if isfield(entity.assets,'centroid_index'),entity.assets=rmfield(entity.assets,'centroid_index');end
    if isfield(entity.assets,'Value_unit'),entity.assets=rmfield(entity.assets,'Value_unit');end
    if isfield(entity.assets,'hazard'),entity.assets=rmfield(entity.assets,'hazard');end
        
    entity.assets = climada_assets_complete(entity.assets);
        
    % create centroid_index already, since we generate the isimip hazard 
    % sets based on this entity and later save time (skip re-encoding)
    entity.assets.centroid_index=1:length(entity.assets.lon);

    if verbose,fprintf('saving entity as %s\n',entity.assets.filename);end
    save(entity.assets.filename,'entity','-v7.3'); % -v7.3 for size...
    
    if check_plot,climada_entity_plot(entity);end
else
    entity=[];
    return
end

end % isimip_ssp2_entity