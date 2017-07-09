function [tc_track,save_file]=isimip_ibtracs_load(basin_name)
% climada isimip ibtracs load tc
% MODULE:
%   isimip
% NAME:
%   isimip_ibtracs_load
% PURPOSE:
%   load isimip ibtracs tropical cyclone (TC) track data
%   Either a single track or all tracks within a folder
%
%   A simple code to load tc_track structure of the repsective basin into
%   memory, see isimip_ibtracs_read for details.
%
%   next call: isimip_tc_hazard_set
% CALLING SEQUENCE:
%   tc_track=isimip_ibtracs_load(basin_name)
% EXAMPLE:
%   tc_track=isimip_ibtracs_load('NA'); % all North Atlantic tracks
%   climada_tc_track_info(tc_track,1); % check plot   
% INPUTS:
%   basin_name: the ocean basin name, such as 'EP', 'NA', 'NI', 'SA', 'SI',
%       'SP', 'WP'. If the .mat file ibtracs.mat does not yet exists,
%       tc_track=isimip_ibtracs_read(basin_name) is called
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   tc_track: a climada TC track structure, see e.g. climada_tc_read_unisys_database
%       plus the fields RadiusMaxWind, EnvironmentalPressure
%   save_file: the file (with path) where the tc_track structure has been
%       saved to and read from
% David N. Bresch, david.bresch@gmail.com, 20161226, intial
%-

save_file='';

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('basin_name','var'),basin_name='';end % OR:

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%

fprintf('DISCONTINUED, use climada_tc_track_load\n');
return

save_file=[climada_global.data_dir filesep 'isimip' filesep ...
    'ibtracs' filesep basin_name filesep 'ibtracs.mat'];

if exist(save_file,'file')
    load(save_file)
else
    [tc_track,save_file]=isimip_ibtracs_read(basin_name,'',1);
end

end % isimip_ibtracs_load
