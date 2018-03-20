function centroids=isimip_admin1_layer(centroids_file,country_ISO3,check_plot,save_it)
% isimip_admin1_layer add admin1 info to isimip centroids
% MODULE:
%   isimip
% NAME:
%   isimip_admin1_layer
% PURPOSE:
%   add admin1 layer to isimip centroids or an entity. The use with an
%   entity might in fact by now be the default use, i.e.
%   entity=isimip_admin1_layer(entity,'ALL')
%
%   search for one admin1, here 'Jura' in Switzerland (CHE):
%   centroids=isimip_admin1_layer('','CHE',2);admin1_name='Jura';
%   for i=1:length(centroids.admin1_list)
%       if strfind(centroids.admin1_list{i,1},admin1_name),pos=i;end;end
%   fprintf('%s has ID %i (at position %i)\n',admin1_name,centroids.admin1_list{pos,2},pos);
%   centroids.admin1_list(pos,:) % to check
%
%   but in case you run a multi-country list, you better use the 3rd entry, i.e.
%   centroids=isimip_admin1_layer('',{'FRA','CHE'},2);admin1_name='Jura | CHE-160';
%   for i=1:length(centroids.admin1_list)
%       if strfind(centroids.admin1_list{i,3},admin1_name),pos=i;end;end
%   fprintf('%s has ID %i (at position %i)\n',admin1_name,centroids.admin1_list{pos,2},pos);
%   centroids.admin1_list(pos,:) % to check, as there is also a 'Jura | FRA-5312' ...
%
%   previous call: <note the most usual previous call here>
%   next call: <note the most usual next function call here>
% CALLING SEQUENCE:
%   res=climada_template(param1,param2);
% EXAMPLE:
%   centroids=isimip_admin1_layer; % TEST, runs it for DEU
%   centroids=isimip_admin1_layer('','DEU',2); % show plot for DEU
%   centroids=isimip_admin1_layer('',{'DEU','FRA'});
%   entity=isimip_admin1_layer(entity,'ALL'); % add admin1 to an entity
% INPUTS:
%   centroids_file: the centroids, if empty (=''), the default isimip centroids are
%       used (i.e. GLB_NatID_grid_0360as_adv_1)
%       if ='ALL', run all countries
%       OR a centroids structure or an entity structure
%       SPECIAL: if centroids_file contains an entity structure for
%           Bangladesh, the special code isimip_BGD_admin is invoked, since
%           Bangladesh added another admin1 region 'Mymensingh'
%   country_ISO3: a single country ISO3 code (like 'USA') or a list of
%       codes (like {'DEU','FRA'}). Default='ALL' (might take time)
% OPTIONAL INPUT PARAMETERS:
%   check_plot: if =1, plot admin1 shapes, default=0
%       if =2, plot all centroids within shapes (might take a lot of time
%       for large countries)
%   save_it: whether we save (=1) or not (=0, default) back to the
%       centroids file
% OUTPUTS:
%   centroids: the centroids, as returned by climada_centroids_load, with
%       fields such as:
%       lon(i), lat(i): longitude and latitude of each centroid i
%       NatID(i): a unique ID for centroid i, see ISO3_list
%       ISO3_list{j,2}: a list with ISO3 country code {:,1} for each NatID in {:,2}
%       and new additional fields:
%       admin1_ID(i): the admin1 index for centroid i, see admin1_name_list
%       admin1_list{j,3}: a list with admin1 name {:,1} for each admin1_ID
%           in {:,2} and admin1 code in {:,3}
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171006, initial
% David N. Bresch, david.bresch@gmail.com, 20171015, small bug fixed (admin1_pos)
% David N. Bresch, david.bresch@gmail.com, 20180306, process entities, too
% David N. Bresch, david.bresch@gmail.com, 20180320, for BGD, call isimip_BGD_admin
%-

%centroids=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('centroids_file','var'),centroids_file='';end
if ~exist('country_ISO3','var'),country_ISO3='';end
if ~exist('check_plot','var'),check_plot=0;end
if ~exist('save_it','var'),save_it=0;end

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
if isempty(centroids_file),centroids_file='GLB_NatID_grid_0360as_adv_1';country_ISO3='DEU';end
if isempty(country_ISO3),country_ISO3{1}='ALL';end
%
module_data_dir=[fileparts(fileparts(which('centroids_generate_hazard_sets'))) filesep 'data'];
admin1_shape_file=[module_data_dir filesep 'ne_10m_admin_1_states_provinces' filesep 'ne_10m_admin_1_states_provinces.shp'];
%
plot_colors={'.r','.g','.b','.k','.m','.y'};

if ~isstruct(centroids_file)
    % load centroids
    centroids=climada_centroids_load(centroids_file);
else
    % contains in fact centroids or an entity, check
    if isfield(centroids_file,'assets')
        % an entity
        
        if isfield(centroids_file.assets,'admin0_ISO3')
            if strcmpi(centroids_file.assets.admin0_ISO3,'BGD') || strcmpi(country_ISO3,'BGD')
                fprintf('SPECIAL for Bangladesh - calling isimip_BGD_admin\n');
                centroids_file=isimip_BGD_admin(centroids_file,'',check_plot);
                centroids=centroids_file;
                return
            end
        end
        
        centroids.lon=centroids_file.assets.lon;
        centroids.lat=centroids_file.assets.lat;
        centroids.filename='ENTITY';
    else
        centroids=centroids_file;
    end
end

% load shape file
admin1_shapes=climada_shaperead(admin1_shape_file);

if ischar(country_ISO3) % convert to cell, if single char
    country_ISO3_tmp=country_ISO3;
    country_ISO3={};country_ISO3{1}=country_ISO3_tmp;
end

if ~strcmpi(country_ISO3{1},'ALL')
    % find all admin1s for all countries as requested
    admin1_shape_i=[]; % init
    for shape_i=1:length(admin1_shapes)
        for country_i=1:length(country_ISO3)
            if strcmp(country_ISO3{country_i},admin1_shapes(shape_i).adm0_a3)
                admin1_shape_i(end+1)=shape_i;
            end
        end % country_i
    end % shape_i
else
    admin1_shape_i=1:length(admin1_shapes);
end

if check_plot,fprintf('NOTE: plotting will take quite some time, be patient\n');end

n_shapes = length(admin1_shape_i);
fprintf('processing %i admin1 shapes\n',n_shapes);
admin1_name_list=cell(1,n_shapes); % init
admin1_name_code_list=cell(1,n_shapes); % init
centroids.admin1_ID=centroids.lon*0; % init
no_centroids_count=0; % init
climada_progress2stdout % init, see terminate below
for admin1_i=1:n_shapes
    
    shape_i=admin1_shape_i(admin1_i);
    admin1_name_list{admin1_i}=admin1_shapes(shape_i).name; % compile list of admin1 names
    admin1_name_code_list{admin1_i}=[admin1_shapes(shape_i).name ...
        ' | ' admin1_shapes(shape_i).adm1_code]; % with code
    
    shape_X=admin1_shapes(shape_i).X;
    shape_Y=admin1_shapes(shape_i).Y;
    if strcmp(admin1_shapes(shape_i).name,'Alaska')
        % SPECIAL case for Alaska (to avoid badly shaped map)
        pos=find(shape_X>100);
        shape_X(pos)=shape_X(pos)-360;
    end
    
    if check_plot % plot admin1 shapes and name them
        plot(shape_X,shape_Y,'-r','LineWidth',1);hold on
        text(admin1_shapes(shape_i).longitude,admin1_shapes(shape_i).latitude,admin1_shapes(shape_i).name);
    end % check_plot
    
    admin1_pos=climada_inpolygon(centroids.lon,centroids.lat,shape_X,shape_Y);
    
    if sum(admin1_pos)>0
        centroids.admin1_ID(admin1_pos)=admin1_i;
        if check_plot>1,plot(centroids.lon(admin1_pos),centroids.lat(admin1_pos),...
                plot_colors{mod(admin1_i,length(plot_colors)-1)+1},'MarkerSize',1);end
    else
        no_centroids_count=no_centroids_count+1;
    end
    
    climada_progress2stdout(admin1_i,n_shapes,10,'admin1 shapes'); % update
    
end % admin1_i
climada_progress2stdout(0) % terminate

if check_plot,set(gcf,'Color',[1 1 1]);end % white figure background
    
if no_centroids_count>0
    fprintf('NOTE: %i admin1 shapes with no centroids\n',no_centroids_count);
end

% reformat
centroids.admin1_list=admin1_name_list';
for i=1:n_shapes
    %centroids.admin1_list{i,1}=strrep(centroids.admin1_list{i,1},'?','_'); % get rid of strange characters
    centroids.admin1_list{i,2}=i;
    centroids.admin1_list{i,3}=admin1_name_code_list{i};
end

centroids.admin1_comment='admin1_name_list(centroids.admin1_ID(i)) is admin1 name for centroids i';

if strcmpi(centroids.filename,'entity')
    % store back to entity
    centroids_file.assets.admin1_ID   = centroids.admin1_ID;
    centroids_file.assets.admin1_list = centroids.admin1_list;
    clear centroids
    centroids=centroids_file; % contains the entity
end

if save_it
    fprintf('saving to %s\n',centroids.filename)
    save(centroids.filename,'centroids',climada_global.save_file_version);
end % save_it

end % isimip_admin1_layer