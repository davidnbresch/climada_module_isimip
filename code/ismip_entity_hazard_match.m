function entity=ismip_entity_hazard_match(entity,hazard,hazard2,check_plot)
% climada template
% MODULE:
%   isimip
% NAME:
%   ismip_entity_hazard_match
% PURPOSE:
%   given an (isimip) entity and a hazard set, join the hazard to the
%   entity, making use of centroid_ID, instead of running a full encoding
%   (which would lead to the same, but take quite long given the 1.6
%   million centroids at 0360as resolution)
%
%   previous call: isimip_gdp_entity, isimip_tc_hazard_set
%   next call: many
% CALLING SEQUENCE:
%   entity=ismip_entity_hazard_match(entity,hazard,hazard2)
% EXAMPLE:
%   entity=ismip_entity_hazard_match('GLB_0360as_entity','Trial1_GB_dkgfdl_20thcal_N_0360as','AUTO',1)
% INPUTS:
%   entity: an isimip entity, either a filename (path automatically
%       figured) or a structure already in memory
%       > promted for if not given
%   hazard: an isimip hazard set, either a filename (path automatically
%       figured) or a structure already in memory
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   hazard2: an second isimip hazard set, either a filename (path automatically
%       figured) or a structure already in memory. This is used if the
%       hazard sets cover northern (N) or southern (S) hemisphere. The code
%       can thus join them.
%       if ='AUTO', the southern hemisphere is inferred automatically (as
%       it is _N_ in the filename is replaced with of _s_). Caution: if the
%       filename contains 'other' _N_, they are replced, too :-)
%       Please note that we take all fields from first hazard, except
%       additional intensities at centroids defined in hazard2.
%   check_plot: =1, plot centroids of assets, hazard and joined hazard
%       =0 no plot (default)
% OUTPUTS:
%   entity: the entity with a field hazard, containing the hazard at
%       exactly the entity centroid's positions
%       There is one new field entity.hazard.stored(centroid_i), which is
%       =1 if we stored intensity for this centroid
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170827, initial
%-

%global climada_global
%if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('entity','var'),entity='';end
if ~exist('hazard','var'),hazard='';end
if ~exist('hazard2','var'),hazard2='';end
if ~exist('check_plot','var'),check_plot=0;end

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%

% load both (does all the handling, i.e. prompting for if not defined...)
entity=climada_entity_load(entity);
if isempty(entity),return;end
hazard=climada_hazard_load(hazard);n_hazards=1;
if isempty(hazard),entity=[];return;end

if strcmpi(hazard2,'AUTO')
    [fP,fN,fE]=fileparts(hazard.filename);
    fN=strrep(fN,'_N_','_S_');
    hazard2=[fP filesep fN fE];
    %     if exist(hazard2,'file')
    %         n_hazards=2;
    %     else
    %         fprintf('WARNING: 2nd hazard %s not found, skipped\n',[fN fE]);
    %     end
end

if ~isstruct(hazard2)
    if exist(hazard2,'file'),n_hazards=2;end
else
    n_hazards=2; % a structure, we assume it's a hazard (checked below)
end

% briefly check consistency

if check_plot
    % to show centroids of both (usually switched off, since boring)
    figure;plot(hazard.lon,hazard.lat,'xy','MarkerSize',3);hold on;
    plot(entity.assets.lon,entity.assets.lat,'.b','MarkerSize',1);
    climada_plot_world_borders;
    title('hazard (yellow), entity assets (blue) and populated (green) centroids')
    xlim([-180 180]);ylim([-90 90])
end

% get rid of possibly exisign hazard
if isfield(entity,'hazard'),entity=rmfield(entity,'hazard');end

n_events=size(hazard.intensity,1);event1=1;
n_centroids=length(entity.assets.lon);
hazard.matrix_density = nnz(hazard.intensity)/numel(hazard.intensity); % recalc

% allocate hazard ( we take all fields from first hazard, except additional
% intensities)
entity.hazard.lon=entity.assets.lon;
entity.hazard.lat=entity.assets.lat;
% NOTE: we keep the entity.hazard.centroid_ID =0 for centroids no
% intensity got stored (see entity.assets.centroid_index for the full list)
entity.hazard.centroid_ID=entity.assets.lon*0; % init
entity.hazard.intensity=spalloc(n_hazards*n_events,n_centroids,ceil(n_hazards*n_events*n_centroids*hazard.matrix_density));
% and all the fields that are define donce for each event
entity.hazard.frequency=hazard.frequency;
entity.hazard.yyyy=hazard.yyyy;
entity.hazard.name=hazard.name;
entity.hazard.comment=hazard.comment;
if isfield(hazard,'windfield_comment'),entity.hazard.windfield_comment=hazard.windfield_comment;end
entity.hazard.peril_ID=hazard.peril_ID;
entity.hazard.units=hazard.units;
entity.hazard.filename=hazard.filename;
entity.hazard.date=hazard.date;
entity.hazard.annotation_str=hazard.annotation_str;
entity.hazard.reference_year=hazard.reference_year;
entity.hazard.orig_years=hazard.orig_years;
entity.hazard.orig_event_count=hazard.orig_event_count;
entity.hazard.orig_event_flag=hazard.orig_event_flag;
entity.hazard.event_count=hazard.event_count;
if isfield(hazard,'ID_no'),entity.hazard.ID_no=hazard.ID_no;end

entity.hazard.stored=entity.hazard.lat*0;

for hazard_i=1:n_hazards;
    fprintf('processing %i events at up to %i centroids\n',size(hazard.intensity,1),n_centroids);
    % allocate check arrays, so we can check correct indexation
    entity.hazard.lon_CHECK=entity.hazard.lon;
    entity.hazard.lat_CHECK=entity.hazard.lat;
    
    % find the intersection of points in entity with points where hazard is defined
    [C,IA,IB] = intersect(entity.assets.centroid_index,hazard.centroid_ID,'stable');
    % returns C such that C=entity.assets.centroid_index(IA)=hazard.centroid_ID(IB)
    if ~isempty(C)
        fprintf('populating %2.0f%% of entity centroids with hazard intensity\n',length(C)/length(entity.assets.lon)*100);
        entity.hazard.intensity(event1:n_events,IA)=hazard.intensity(:,IB); % (event,centroid)
        entity.hazard.centroid_ID(IA)=hazard.centroid_ID(IB); % for records
        entity.hazard.stored(IA)=1; % to note the centroids we stored intensity for

        % check (uisng lat/lon and making sure the match for all the
        % centroids we stored intensities for)
        entity.hazard.lon_CHECK(IA)=hazard.lon(IB);
        entity.hazard.lat_CHECK(IA)=hazard.lat(IB);
        lola_check=abs(sum(entity.hazard.lon-entity.hazard.lon_CHECK))+abs(sum(entity.hazard.lat-entity.hazard.lat_CHECK));
        if lola_check>10*eps,fprintf('WARNING: matching might be wrong (lola checksum %f\n',lola_check);end
    else
        fprintf('WARNING: no hazard at entity centroids\n');
    end
        
    if n_hazards==2 && hazard_i==1 % load 2nd hazard 
        event1=n_events+1;
        hazard=climada_hazard_load(hazard2);
        if isempty(hazard),break;end
        n_events=n_events+size(hazard.intensity,1);
        fprintf('> dealing with 2nd hazard %s (total %i events)\n',[fN fE],n_events);
        entity.hazard.comment=[entity.hazard.comment ' + ' hazard.comment];
        if isfield(hazard,'windfield_comment'),entity.hazard.windfield_comment=[entity.hazard.windfield_comment ' + ' hazard.windfield_comment];end
        if check_plot,plot(hazard.lon,hazard.lat,'oy','MarkerSize',3);end % add centroids of 2nd hazard
    end
    
end % hazard_i

entity.hazard.intensity=entity.hazard.intensity(1:n_events,:); % fix exact size

% complete hazard
entity.hazard.matrix_density = nnz(entity.hazard.intensity)/numel(entity.hazard.intensity); % recalc
%hazard.fraction=spones(hazard.intensity); % fraction 100% WE DO NOT ADD for SAVINF MEMORY

if check_plot % to show centroids populated with hazard
    stored_pos=find(entity.hazard.stored);
    plot(entity.hazard.lon(stored_pos),entity.hazard.lat(stored_pos),'.g','MarkerSize',1);
    climada_plot_world_borders;
    title('hazard (yellow), entity assets (blue) and populated (green) centroids')
    xlim([-180 180]);ylim([-90 90])
    % saveas(gcf,'assets_hazard_centroids.png','png')
end

end % ismip_entity_hazard_match