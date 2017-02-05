function [YDS,EDS_stats]=isimip_YDS_calc(entity,hazard,params)
% climada isimip YDS EDS calc
% MODULE:
%   isimip
% NAME:
%   isimip_YDS_calc
% PURPOSE:
%   calculate year damage set (YDS) for isimip assets and hazards
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
%   EDS=isimip_YDS_calc(entity,hazard)
% EXAMPLE:
%   one country:
%   entity=isimip_gdp_entity('DEU') % create single country entity
%   entity=climada_entity_load('DEU_entity') % load DEU entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create DEU FL hazard set
%   [YDS,EDS_stats]=isimip_YDS_calc(entity,hazard,1);
%
%   entity=isimip_gdp_entity('USA') % create single country entity
%   entity=climada_entity_load('USA_entity') % load entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create USA FL hazard set
%   [YDS,EDS_stats]=isimip_YDS_calc(entity,hazard);
%
%   list of countries:
%   entity=isimip_gdp_entity({'DEU','ITA','FRA'}) % create entity for DEU ITA and FRA
%   entity=climada_entity_load('DEUITAFRA_entity') % load DEU entity
%   hazard=isimip_flood_load('global_FL.nc','auto',entity,0); % create DEU FL hazard set
%   [YDS,EDS_stats]=isimip_YDS_calc(entity,hazard);
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
%       separately. Default=1 (plot), no plot if =0
%    markersize: the marker size for the check plot, default=3
% OUTPUTS:
%   YDS: the year damage set (YDS), see e.g. climada_EDS2YDS for details
%   PLUS the fields
%       YDS.lon and YDS.at lat/lon for ED_at_centroid (if YDS has more than
%           one element)
%       NOTE: to plot the expected damage at each centroid for YDS i, use e.g.
%       a.assets.lon=YDS(i).lon;a.assets.lat=YDS(i).lat;
%       a.assets.Value=YDS(i).ED_at_centroid';climada_entity_plot(a);
%       title(YDS(i).comment)
%   EDS_stats: some statistics (i.e. entity and hazard year matches etc.)
%       entity_yyyy: the years as in the entity (all)
%       hazard_yyyy: the years there has been a matching hazard (=NaN no
%           hazard for a particular corresponding year in entity_yyyy)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161120, initial
% David N. Bresch, david.bresch@gmail.com, 20170107, damage_at_centroid
%-

YDS=[];EDS_stats=[]; % init output

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

% PARAMETERS
%
% set default values (see header for details)
if isempty(params.check_plot),params.check_plot=1;end
if isempty(params.markersize),params.markersize=3;end
%
silent_mode=0;

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

EDS_stats.entity_yyyy=entity.assets.Values_yyyy;
EDS_stats.hazard_yyyy=EDS_stats.entity_yyyy+NaN; % init

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

% figure the matching years
valid_time_i=[]; % init
for time_i=1:n_times
    hazard_yyyy_pos=find(hazard.yyyy==entity.assets.Values_yyyy(time_i), 1);
    if ~isempty(hazard_yyyy_pos)
        valid_time_i=[valid_time_i time_i];
    end
end % time_i

n_times=length(valid_time_i);

% pre-fill the YDS
for country_i=1:n_countries
    YDS(country_i).reference_year=-999; % init
    YDS(country_i).Value=0; % init
    YDS(country_i).ED=0; % init
    [~,list_pos]=intersect(cell2mat(entity.assets.ISO3_list(:,2)),unique_NatId(country_i));
    YDS(country_i).comment=entity.assets.ISO3_list{list_pos,1};
    YDS(country_i).assets.filename=entity.assets.filename; % as expected by eg climada_EDS_DFC
    country_info(country_i).pos=find(entity.assets.NatId==unique_NatId(country_i));
    YDS(country_i).lon=entity.assets.lon(country_info(country_i).pos);
    YDS(country_i).lat=entity.assets.lat(country_info(country_i).pos);
    YDS(country_i).damage_at_centroid=zeros(n_times,length(YDS(country_i).lon)); % init
end % country_i
next_YDS_i=1; % init

% template for-loop with waitbar or progress to stdout
t0       = clock;
mod_step = 1; % first time estimate after 10 events, then every 100 (see below)
format_str='%s';

msgstr   = sprintf('processing %i years (%i countries)',n_times,n_countries);
fprintf('%s\n',msgstr);

for time_i=1:n_times
    
    time_ii=valid_time_i(time_i);
    
    % your calculations here
    entity.assets.Value=entity.assets.Values(time_ii,:); % the assets for time i
    
    % find corresponding hazard
    hazard_yyyy_pos=find(hazard.yyyy==entity.assets.Values_yyyy(time_ii));
    
    if ~isempty(hazard_yyyy_pos)
        EDS_stats.hazard_yyyy(time_ii)=EDS_stats.entity_yyyy(time_ii);
        
        % copy only the events for the corresponding year
        temp_hazard.intensity=hazard.intensity(hazard_yyyy_pos,:);
        temp_hazard.fraction =hazard.fraction(hazard_yyyy_pos,:);
        temp_hazard.frequency=hazard.frequency(hazard_yyyy_pos)*0+1; % all one (i.e. we sum up)
        
        EDS=climada_EDS_calc(entity,temp_hazard,'',0,2);
        
        for country_i=1:n_countries
            CP=country_info(country_i).pos; % just for shorter code
            
            YDS(country_i).damage(next_YDS_i)=sum(EDS.ED_at_centroid(CP));
            YDS(country_i).damage_at_centroid(next_YDS_i,:)=EDS.ED_at_centroid(CP)';
            YDS(country_i).Value=YDS(country_i).Value+entity.assets.Value(CP);
            
            YDS(country_i).event_ID(next_YDS_i)=time_ii;
            YDS(country_i).frequency(next_YDS_i)=1;
            YDS(country_i).yyyy(next_YDS_i)=EDS_stats.entity_yyyy(time_ii);
            
            YDS(country_i).peril_ID=EDS.peril_ID;
            YDS(country_i).hazard=EDS.hazard;
            YDS(country_i).Value_unit=EDS.Value_unit;
            YDS(country_i).damagefunctions=EDS.damagefunctions;
            YDS(country_i).annotation_name=EDS.annotation_name;
        end % country_i
        
        next_YDS_i=next_YDS_i+1;
    end % ~isempty(hazard_yyyy_pos)
    
    % the progress management
    if mod(time_i,mod_step)==0
        mod_step          = 2;
        t_elapsed_time   = etime(clock,t0)/time_i;
        times_remaining  = n_times-time_i;
        t_projected_sec   = t_elapsed_time*times_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i years)',t_projected_sec,   time_i,n_times);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i years)',t_projected_sec/60,time_i,n_times);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end % time_i
fprintf(format_str,''); % move carriage to begin of line

% finish YDS
n_years=length(YDS(1).damage);
for country_i=1:n_countries
    YDS(country_i).reference_year=max(min(EDS_stats.entity_yyyy),min(EDS_stats.hazard_yyyy));
    YDS(country_i).orig_event_flag= YDS(country_i).damage*0+1;
    YDS(country_i).orig_year_flag = YDS(country_i).damage*0+1;
    YDS(country_i).Value          = YDS(country_i).Value/n_years; % take average
    YDS(country_i).ED             = sum(YDS(country_i).damage)/n_years; % take average
    YDS(country_i).frequency      = YDS(country_i).frequency/n_years; % per year
end % country_i

if params.check_plot
    m=ceil(sqrt(length(YDS)));n=m;
    fprintf('plotting');
    for i=1:length(YDS)
        subplot(m,n,i)
        damage_at_centroid=sum(YDS(i).damage_at_centroid,1)/n_years;
        if sum(damage_at_centroid)>0 % something to plot
            fprintf(' %s',YDS(i).comment);
            a.assets.lon=YDS(i).lon;a.assets.lat=YDS(i).lat;
            a.assets.Value=damage_at_centroid;
            climada_entity_plot(a,params.markersize);hold on;title([YDS(i).comment ' ' YDS(i).peril_ID ' ex. damage']);
        else
            title([YDS(i).comment ' ' YDS(i).peril_ID ' damage ZERO']);
        end
    end % country_i
    fprintf('\n');
end % params.check_plot

end % isimip_YDS_calc