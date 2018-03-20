function entity=isimip_BGD_admin(entity,shapefile_name,check_plot)
% climada template
% MODULE:
%   isimip
% NAME:
%   isimip_BGD_admin
% PURPOSE:
%   add admin 1 layer to an entity (for Bangladesh at the moment),
%   see shapefile_name in INPUTS and in PARAMETERS
%
%   shape file from:
%   https://data.humdata.org/dataset/administrative-boundaries-of-bangladesh-as-of-2015
%
%   previous call: isimip_gdp_entity
%   next call: many
% CALLING SEQUENCE:
%   entity=isimip_BGD_admin(entity,shapefile_name,check_plot)
% EXAMPLE:
%   entity=isimip_BGD_admin(isimip_gdp_entity('BGD','0150as'))
% INPUTS:
%   entity:  a climada entity, see e.g. climada_entity_load
%       > prompted for if not given
% OPTIONAL INPUT PARAMETERS:
%   shapefile_name: the name of a shape file (.shp) with admin information,
%       see PARAMETRES for default (currently for Bangladesh)
%       fields used are X, Y and adm1_en
%   check_plot: if=1, plot shapes, =0 not (default)
%       if >1, plot also centroids within admin1 boundaries
% OUTPUTS:
%   entity: the entity structure as on input, just with fields added:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20180320, initial
%-

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments and to set default values
if ~exist('entity','var'),entity=[];end
if ~exist('shapefile_name','var'),shapefile_name='';end
if ~exist('check_plot','var'),check_plot=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% set default shape file name
if isempty(shapefile_name),shapefile_name=[module_data_dir filesep 'centroids' filesep 'bgd_bnda_adm1_2015/bgd_bnda_adm1_2015.shp'];end
%shapefile_name='/Users/bresch/Documents/_GIT/climada_applications_div/isimip_BGD/admin_layers/bgd_bnda_adm1_2015/bgd_bnda_adm1_2015.shp';
%
plot_colors={'.r','.g','.b','.k','.m','.y'};

entity=climada_entity_load(entity);

if exist(shapefile_name,'file')
    shapes=climada_shaperead(shapefile_name);
else
    fprintf('ERROR: shapefile %s not found\n',shapefile_name);
end

% template for-loop with progress to stdout
n_shapes = length(shapes);
fprintf('processing %i admin shapes for %i centroids\n',n_shapes,length(entity.assets.lon));
no_centroids_count=0;
entity.assets.admin1_ID=entity.assets.lon*0;
admin1_name_list=cell(1,n_shapes); % init
admin1_name_code_list=cell(1,n_shapes); % init

climada_progress2stdout(-1,[],1) % init with mod_step 1, see terminate below
for shape_i=1:n_shapes
    
    admin1_name_list{shape_i}=shapes(shape_i).adm1_en; % compile list of admin1 names
    admin1_name_code_list{shape_i}=[admin1_name_list{shape_i} ...
        ' | ' num2str(shape_i)]; % with code
    
    shape_X=shapes(shape_i).X;
    shape_Y=shapes(shape_i).Y;
    
    if check_plot % plot admin1 shapes and name them
        plot(shape_X,shape_Y,'-r','LineWidth',1);hold on
    end % check_plot
    
    admin1_pos=climada_inpolygon(entity.assets.lon,entity.assets.lat,shape_X,shape_Y);
    
    if sum(admin1_pos)>0
        entity.assets.admin1_ID(admin1_pos)=shape_i; % simply the position
        if check_plot>1,plot(entity.assets.lon(admin1_pos),entity.assets.lat(admin1_pos),...
                plot_colors{mod(shape_i,length(plot_colors)-1)+1},'MarkerSize',1);end
    else
        no_centroids_count=no_centroids_count+1;
    end
    
    climada_progress2stdout(shape_i,n_shapes,1,'shapes'); % update
end % shape_i
climada_progress2stdout(0) % terminate

if check_plot % plot admin1 shapes and name them
    axis equal
    xlim0=xlim;ylim0=ylim;
    climada_plot_world_borders
    xlim(xlim0);ylim(ylim0);
    for shape_i=1:n_shapes
        shape_X=shapes(shape_i).X;shape_Y=shapes(shape_i).Y;
        text(mean(shape_X(~isnan(shape_X))),mean(shape_Y(~isnan(shape_Y))),shapes(shape_i).adm1_en);
    end
end % check_plot

% reformat
entity.assets.admin1_list=admin1_name_list';
for i=1:n_shapes
    entity.assets.admin1_list{i,2}=i;
    entity.assets.admin1_list{i,3}=admin1_name_code_list{i};
end

end % isimip_BGD_admin