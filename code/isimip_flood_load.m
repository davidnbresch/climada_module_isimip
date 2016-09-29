function hazard=isimip_flood_load(flood_filename,entity,check_plot)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_flood_load
% PURPOSE:
%   load .mat file with flood footprints
%
%   next call: 
% CALLING SEQUENCE:
%   hazard=isimip_flood_load(flood_filename,entity,check_plot)
% EXAMPLE:
%   entity=climada_entity_load('USA_UnitedStates_Florida');
%   hazard=isimip_flood_load('fldfrc_max.nc',entity,1);
% INPUTS:
%   flood_filename: filename of the .nc file with the flood footprints
%       tracks, default folder is ..climada_data/isimip
%       > promted for if not given
%   entity: an entity struct to interpolate the flood footprints to, see
%       climada_entity_load and climada_entity_read for a description
% OPTIONAL INPUT PARAMETERS:
%   check_plot: whether show a check plot (=1, default), or not (=0)
%       Note that plotting might often take longer than the full
%       conversion...
% OUTPUTS:
%   hazard: a climada hazard structure, see manual
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160929, initial
%-

hazard=[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('flood_filename','var'),flood_filename='';end
if ~exist('entity','var'),entity='';end
if ~exist('check_plot','var'),check_plot=1;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];

flood_fraction_filename='fldfrc_max.nc';
flood_depth_filename='flddph_max.nc';

% template to prompt for track_filename if not given
if isempty(flood_filename) % local GUI
    flood_filename=[isimip_data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(flood_filename, 'Open:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        flood_filename=fullfile(pathname,filename);
    end
end

entity=climada_entity_load(entity);

% complete path, if missing
[fP,fN,fE]=fileparts(flood_fraction_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_fraction_filename=[fP filesep fN fE];
if ~exist(flood_fraction_filename,'file')
    fprintf('Error: %s nt found\n',flood_fraction_filename);
end

% complete path, if missing
[fP,fN,fE]=fileparts(flood_depth_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_depth_filename=[fP filesep fN fE];
if ~exist(flood_depth_filename,'file')
    fprintf('Error: %s nt found\n',flood_depth_filename);
end

fprintf('reading from %s and %s ...',flood_fraction_filename,flood_depth_filename);
nc.lon      = ncread(flood_fraction_filename,'lon');
nc.lat      = ncread(flood_fraction_filename,'lat');
nc.time     = ncread(flood_fraction_filename,'time');
dlat=1586-1336;
dlon=2469-2205;
% if troubles, use ncinfo(flood_fraction_filename,'var')
nc.fraction = ncread(flood_fraction_filename,'var',[2205 1336 1],[dlon dlat 1]); % time, lat, lon
nc.depth    = ncread(flood_depth_filename,   'var',[2205 1336 1],[dlon dlat 1]); % time, lat, lon

% fraction 0..1, depth in meters

fprintf(' done\n');

% interpolate to entity, then call sven.willner@pik-potsdam.de 


end % isimip_flood_load


if check_plot
    fprintf('plotting %i tracks ...',length(tc_track))
   %
    climada_plot_world_borders
    fprintf('done\n')
end % check_plot

fprintf('saving hazard as %s\n',save_filename)
save(save_filename,'tc_track');

end % isimip_load_tc_tracks