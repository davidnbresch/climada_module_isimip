function [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard,params)
% climada isimip YDS EDS calc
% MODULE:
%   isimip
% NAME:
%   isimip_YDS_calc
% PURPOSE:
%   calculate year damage set (YDS) for isimip assets and hazards
%
%   TODO: matching_Natcat-damages_ibtracs_1980-2014
%
%   for each year where there are assets, find the corresponding events and
%   calculate the event damage set (EDS) for this particular year, then
%   store the year damage set (YDS, i.e. summed over events within given year).
%
%   If the entity does contain more than one (isimip) country, create one
%   YDS for each country separately (based upon entity.assets.NatId).
%
%   previous call: isimip_gdp_entity and e.g. isimip_flood_load
%   next call: isimip_YDS_table
% CALLING SEQUENCE:
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard,params)
% EXAMPLE:
%   one country:
%   entity=isimip_gdp_entity('DEU') % create single country entity
%   entity=climada_entity_load('DEU_entity') % load DEU entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create DEU FL hazard set
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard,1);
%
%   entity=isimip_gdp_entity('USA') % create single country entity
%   entity=climada_entity_load('USA_entity') % load entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create USA FL hazard set
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard);
%
%   entity=isimip_gdp_entity('USA') % create single country entity
%   entity=climada_entity_load('USA_entity') % load entity
%   hazard=climada_hazard_load('GLB_glb_TC'); % see isimip_tc_hazard_set
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard);
%
%   entity=isimip_gdp_entity({'PRI','DOM','BRB','CUB'}) % create single country entity
%   entity=climada_entity_load('BRBCUBDOMPRI') % load entity
%   hazard=climada_hazard_load('GLB_0360as_TC'); % see isimip_tc_hazard_set
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard);
%
%   list of countries:
%   entity=isimip_gdp_entity({'DEU','ITA','FRA'}) % create entity for DEU ITA and FRA
%   entity=climada_entity_load('DEUITAFRA_entity') % load DEU entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create DEU FL hazard set
%   [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard);
% INPUTS:
%   entity: an isimip entity (i.e. an entity with entity.assets.Values for
%       many years). Needs to contain the following additional a fields:
%       (a so-called isimip entity has that, e.g. if created by isimip_gdp_entity) 
%       entity.assets.NatId(centroid_i): to link each centroid to one country  
%       entity.assets.Values(year_i,centroid_i): the assets for year_i
%       entity.assets.Value_yyyy(year_i): the year assets are valid for
%       If ='params', return all default parameters in YDS, params=isimip_YDS_calc('params')
%   hazard: a hazard event set (either a regular one or an isimip special
%       one)
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%    check_plot: if=1, plot exected damage at each centroid for each country
%       separately. Default=0 (no plot)
%    markersize: the marker size for the check plot, default=3
%    damage_function_regions_file: the filename (with path) where the damage
%       functions regions (groups of countries) are defined. A .csv file
%       with columns ISO (3-char), ID (int), Reg_ID (int) and Reg_name (3-char)
% OUTPUTS:
%   YDS: the year damage set (YDS), see e.g. climada_EDS2YDS for details
%    PLUS the fields
%       YDS.lon and YDS.lat lat/lon for ED_at_centroid (if YDS has more than
%           one element)
%       NOTE: to plot the expected damage at each centroid for YDS i, use e.g.
%       a.assets.lon=YDS(i).lon;a.assets.lat=YDS(i).lat;
%       a.assets.Value=YDS(i).ED_at_centroid';climada_entity_plot(a);
%       title(YDS(i).comment)
%   EDS(country_i): the event damage set (EDS) for country i
%   stats: some statistics (i.e. entity and hazard year matches etc.)
%       entity_yyyy: the years as in the entity (all)
%       hazard_yyyy: the years there has been a matching hazard (=NaN no
%           hazard for a particular corresponding year in entity_yyyy)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161120, initial
% David N. Bresch, david.bresch@gmail.com, 20170107, damage_at_centroid
% David N. Bresch, david.bresch@gmail.com, 20170216, EDS(country_i) returned
% David N. Bresch, david.bresch@gmail.com, 20170217, simplified
%-

YDS=[];EDS=[];stats=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('entity','var'),         entity=[];end
if ~exist('hazard','var'),         hazard=[];end
if ~exist('params','var'),         params=struct;end

% check for some parameter fields we need
if ~isfield(params,'check_plot'),params.check_plot=[];end
if ~isfield(params,'markersize'),params.markersize=[];end
if ~isfield(params,'damage_function_regions_file'),params.damage_function_regions_file=[];end

% PARAMETERS
%
% set default values (see header for details)
if isempty(params.check_plot),params.check_plot=0;end
if isempty(params.markersize),params.markersize=3;end
if isempty(params.damage_function_regions_file),...
        params.damage_function_regions_file=[climada_global.data_dir filesep 'isimip' filesep 'TC-damage-function-regions.csv'];end

%
silent_mode=0;
%
climada_global.damage_at_centroid=1;

if strcmpi(entity,'params'),YDS=params;return;end % special case, return the full params strcture

% check/process input
entity = climada_entity_load(entity); % prompt for entity if not given
hazard = climada_hazard_load(hazard); % prompt for hazard if not given
if isempty(hazard) || isempty(entity),return;end
hazard = climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

if ~isfield(entity.assets,'NatId')
    fprintf('Warning: no field entity.assets.NatId, might not be an isimip entity\n'); 
end

% figure number of years to process
if ~isfield(entity.assets,'Values_yyyy')
    fprintf('ERROR: no field entity.assets.Values_yyyy, not an isimip entity, aborted\n');
    return
end

n_times=length(entity.assets.Values_yyyy);

% encode assets of entity once more, just to be sure
if ~isfield(entity.assets,'centroid_index')
    if ~silent_mode,fprintf('Encoding entity assets to hazard ... ');end
    entity = climada_assets_encode(entity,hazard);
elseif ~all(diff(entity.assets.centroid_index) == 1) && climada_global.re_check_encoding
    if ~silent_mode,fprintf('Encode entity assets once more ...');end
    entity = climada_assets_encode(entity,hazard);
end

% make a local copy of all the small fields
temp_hazard.peril_ID        =hazard.peril_ID;
temp_hazard.units           =hazard.units;
temp_hazard.reference_year  =hazard.reference_year;
temp_hazard.lon             =hazard.lon;
temp_hazard.lat             =hazard.lat;
temp_hazard.centroid_ID     =hazard.centroid_ID;
temp_hazard.orig_years      =hazard.orig_years;
temp_hazard.orig_event_count=hazard.orig_event_count;
temp_hazard.event_count     =hazard.event_count;
temp_hazard.event_ID        =hazard.event_ID;
temp_hazard.orig_event_flag =hazard.orig_event_flag;
temp_hazard.yyyy            =hazard.yyyy;
temp_hazard.filename        =hazard.filename;
temp_hazard.date            =hazard.date;
temp_hazard.comment         =hazard.comment;
temp_hazard.filename        =hazard.filename;
        
% figure the number of countries in the entity
try
    unique_NatId=unique(entity.assets.NatId);
    n_countries=length(unique_NatId);
catch
    fprintf('Warning: not able to infer number of countries, assuming only one\n');
    n_countries=1; % default to one
end

% just some info we store
stats.entity_yyyy=entity.assets.Values_yyyy;
stats.hazard_yyyy=stats.entity_yyyy+NaN; % init
stats.ISO3_list=entity.assets.ISO3_list;

% figure the matching years
valid_time_i=[]; % init
for time_i=1:n_times
    hazard_yyyy_pos=find(hazard.yyyy==entity.assets.Values_yyyy(time_i), 1);
    if ~isempty(hazard_yyyy_pos)
        valid_time_i=[valid_time_i time_i];
    end
end % time_i

n_times=length(valid_time_i);

for country_i=1:n_countries
    % pre-fill the YDS
    YDS(country_i).reference_year     = entity.assets.Values_yyyy(valid_time_i(1)); % first year
    YDS(country_i).Value              = 0; % init
    YDS(country_i).ED                 = 0; % init
    [~,list_pos] = intersect(cell2mat(entity.assets.ISO3_list(:,2)),unique_NatId(country_i));
    YDS(country_i).comment            = entity.assets.ISO3_list{list_pos,1};
    YDS(country_i).assets.filename    = entity.assets.filename; % as expected by eg climada_EDS_DFC
    country_info(country_i).pos = find(entity.assets.NatId==unique_NatId(country_i));
    YDS(country_i).lon                = entity.assets.lon(country_info(country_i).pos);
    YDS(country_i).lat                = entity.assets.lat(country_info(country_i).pos);
    YDS(country_i).damage_at_centroid = zeros(n_times,length(YDS(country_i).lon)); % init
   
    % pre-fill the EDS
    EDS(country_i).reference_year     = YDS(country_i).reference_year;
    EDS(country_i).event_ID           = hazard.event_ID;
    EDS(country_i).frequency          = hazard.frequency;
    EDS(country_i).orig_event_flag    = hazard.orig_event_flag;
    EDS(country_i).peril_ID           = hazard.peril_ID;
    EDS(country_i).damage             = hazard.event_ID*0; % init
    EDS(country_i).ED_at_centroid     = zeros(1,length(YDS(country_i).lon)); % init
    EDS(country_i).damage_at_centroid = []; % init
    EDS(country_i).Value              = 0; % init
    EDS(country_i).ED                 = 0; % init
    EDS(country_i).comment            = entity.assets.ISO3_list{list_pos,1};
    EDS(country_i).annotation_name    = EDS(country_i).comment;
    EDS(country_i).assets.lon         = YDS(country_i).lon; % init
    EDS(country_i).assets.lat         = YDS(country_i).lat; % init
    EDS(country_i).assets.filename    = EDS(country_i).comment;
    
end % country_i
next_YDS_i=1; % init

fprintf('processing %i years (%i countries)\n',n_times,n_countries);
climada_progress2stdout    % init, see terminate below

for time_i=1:n_times
    
    time_ii=valid_time_i(time_i);
    
    % your calculations here
    entity.assets.Value=entity.assets.Values(time_ii,:); % the assets for time i
    
    % find corresponding hazard
    hazard_yyyy_pos=find(hazard.yyyy==entity.assets.Values_yyyy(time_ii));
    
    if ~isempty(hazard_yyyy_pos)
        stats.hazard_yyyy(time_ii)=stats.entity_yyyy(time_ii);
        
        % copy only the events for the corresponding year
        temp_hazard.intensity       =hazard.intensity(hazard_yyyy_pos,:);
        temp_hazard.fraction        =hazard.fraction(hazard_yyyy_pos,:);
        temp_hazard.frequency       =hazard.frequency(hazard_yyyy_pos)*0+1; % all one (i.e. we sum up)
        
        temp_hazard.orig_event_count=length(hazard_yyyy_pos);
        temp_hazard.event_count     =length(hazard_yyyy_pos);
        temp_hazard.event_ID        =hazard.event_ID(hazard_yyyy_pos);
        temp_hazard.orig_event_flag =hazard.orig_event_flag(hazard_yyyy_pos);
        temp_hazard.yyyy            =hazard.yyyy(hazard_yyyy_pos);
        temp_hazard.orig_years      =1;
                
        EDS_temp=climada_EDS_calc(entity,temp_hazard,'',0,2);
                
        for country_i=1:n_countries
            CP=country_info(country_i).pos; % just for shorter code
            
            YDS(country_i).damage(next_YDS_i)=sum(EDS_temp.ED_at_centroid(CP));
            YDS(country_i).damage_at_centroid(next_YDS_i,:)=EDS_temp.ED_at_centroid(CP)';
            YDS(country_i).Value=YDS(country_i).Value+entity.assets.Value(CP);

            YDS(country_i).event_ID(next_YDS_i)=time_ii;
            YDS(country_i).frequency(next_YDS_i)=1;
            YDS(country_i).yyyy(next_YDS_i)=stats.entity_yyyy(time_ii);
            
            YDS(country_i).peril_ID=EDS_temp.peril_ID;
            YDS(country_i).hazard=EDS_temp.hazard;
            YDS(country_i).Value_unit=EDS_temp.Value_unit;
            YDS(country_i).damagefunctions=EDS_temp.damagefunctions;
            YDS(country_i).annotation_name=EDS_temp.annotation_name;
            
            EDS(country_i).hazard=EDS_temp.hazard;
            EDS(country_i).Value_unit=EDS_temp.Value_unit;
            EDS(country_i).damagefunctions=EDS_temp.damagefunctions;

            EDS(country_i).damage(hazard_yyyy_pos)=sum(EDS_temp.damage_at_centroid(CP,:)); % damage_at_centroid(centroid_i,event_i)
            EDS(country_i).ED_at_centroid=EDS(country_i).ED_at_centroid+EDS_temp.ED_at_centroid(CP)'; % ED_at_centroid(centroid_i)
            EDS(country_i).damage_at_centroid(CP,hazard_yyyy_pos)=EDS_temp.damage_at_centroid(CP,:); % damage_at_centroid(centroid_i,event_i)
            EDS(country_i).Value=EDS(country_i).Value+sum(entity.assets.Value(CP));
            
        end % country_i
        
        next_YDS_i=next_YDS_i+1;
    end % ~isempty(hazard_yyyy_pos)
    
    climada_progress2stdout(time_i,n_times,1,'years'); % the progress management
    
end % time_i
climada_progress2stdout(0) % terminate

n_years=length(YDS(1).damage);
for country_i=1:n_countries
    
    % finish YDS
    YDS(country_i).reference_year = max(min(stats.entity_yyyy),min(stats.hazard_yyyy));
    YDS(country_i).orig_event_flag= YDS(country_i).damage*0+1;
    YDS(country_i).orig_year_flag = YDS(country_i).damage*0+1;
    YDS(country_i).Value          = YDS(country_i).Value/n_years; % take average
    YDS(country_i).ED             = sum(YDS(country_i).damage)/n_years; % take average
    YDS(country_i).frequency      = YDS(country_i).frequency/n_years; % per year
    
    % finish EDS
    EDS(country_i).Value          = EDS(country_i).Value/n_years; % average Value over years
    EDS(country_i).ED             = sum(EDS(country_i).damage)/n_years; % take average

end % country_i

if params.check_plot
    n_countries=length(YDS);
    m=ceil(sqrt(n_countries));
    n=max(m-1,1);if m*n<n_countries,n=m;end
    fprintf('plotting');
    for i=1:n_countries
        subplot(m,n,i)
        damage_at_centroid=sum(YDS(i).damage_at_centroid,1)/n_years;
        if sum(damage_at_centroid)>0 % something to plot
            fprintf(' %s',YDS(i).comment);
            a.assets.lon=YDS(i).lon;a.assets.lat=YDS(i).lat;
            a.assets.Value=damage_at_centroid;a.assets.filename=YDS(i).comment;
            entity_plot_params.cbar_ylabel='expected damage';
            climada_entity_plot(a,params.markersize,entity_plot_params);hold on;title([YDS(i).comment ' ' YDS(i).peril_ID]);
        else
            title([YDS(i).comment ' ' YDS(i).peril_ID ' damage ZERO']);
        end
    end % country_i
    fprintf('\n');
end % params.check_plot

damage_function_regions=climada_csvread(params.damage_function_regions_file);
% ISO: {1x230 cell}
% ID: [1x230 double]
% Reg_ID: [1x230 double]
% Reg_name: {1x230 cell}

params.damage_data_file=[climada_global.data_dir filesep 'isimip' filesep 'matching_Natcat-damages_ibtracs_1980-2014.csv'];
% commas within fields
damage_data=climada_csvread(params.damage_data_file);

params.price_deflator_file=[climada_global.data_dir filesep 'isimip' filesep 'GDP_deflator_converted_base2005_1969-2016_source_BEA.csv'];
deflator_data=climada_csvread(params.price_deflator_file);
deflator_value=1./(deflator_data.GDP_deflator_base2005/100);

end % isimip_YDS_calc