function [flddph_filename,fldfrc_filename,fld_path]=isimip_get_flood_filename(isimip_simround, ghm, forcing, protection, time_period)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_get_flood_filename
% PURPOSE:
%   provide a unique and consistent name for the flood file (.nc) from
%   isimip CaMa-Flood simulations. This avoid using different conventions
%   or repeating code too much.
%
%   next call: isimip_flood_load
% CALLING SEQUENCE:
%   hazard_filename=isimip_get_flood_filename(flood_filename,entity,check_plot,isimip_data_subdir,years_range)
% EXAMPLE:
%   isimip_simround='2b';
%   flood_filename='merged_LPJmL_miroc5_historical_flopros_gev_0.1.nc';
%   entity=climada_entity_load('USA_UnitedStates_Florida');
%   hazard_filename=isimip_get_flood_hazard_filename(flood_filename,entity,isimip_data_subdir)
% INPUTS:
%   flood_filename: filename of the .nc file with the flood
%       footprints, default folder is ..climada_data/isimip/FL
%       > promted for if not given
%       fraction (variable name 'fldfrc') is in the range 0..1
%       depth (variable name 'flddph') in units of meters [m]
%       there should be one event per year (i.e., yearly maxima)
%   entity: an entity struct to interpolate the flood footprints to, see
%       climada_entity_load and climada_entity_read for a description
% OPTIONAL INPUT PARAMETERS:
%   isimip_data_subdir: the sub-directory within the isimip folder within
%       climada_data/isimip where the raw isimip data (NetCDF file
%       flood_filename) is located. This does not affect where the hazard is
%       saved. If not specified, the file is assumed to be located directly
%       within the climada_data/isimip folder.
%   years_range: vector of length 2 containing the first and the last year
%       to be loaded from the netcdf file. If empty or [0 0], loads all data.
% OUTPUTS:
%   hazard_filename: a file name for the .mat hazard file on entity grid.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180118, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('isimip_simround','var'),        isimip_simround=         '';end
if ~exist('ghm','var'),                    ghm=                     '';end
if ~exist('forcing','var'),                forcing=                 '';end
if ~exist('protection','var')              protection=             '0';end
if ~exist('time_period','var'),            time_period=   'historical';end

isimip_data_dir = [climada_global.data_dir filesep 'isimip'];

% define files and path
fld_path=[isimip_data_dir filesep isimip_simround];
if isimip_simround=='2a'
    flddph_filename=['flddph_' ghm '_' forcing '_' protection '_gev_0.1.nc'];
    fldfrc_filename=['fldfrc_' ghm '_' forcing '_' protection '_gev_0.1.nc'];
elseif isimip_simround=='2b'
    flddph_filename=['flddph_' ghm '_' forcing '_' time_period '_' protection '_gev_0.1.nc'];
    fldfrc_filename=['fldfrc_' ghm '_' forcing '_' time_period '_' protection '_gev_0.1.nc'];
elseif isimip_simround==''
    fprintf('No value in isimip_simround, filename will be unknown');
else
    fprintf('Unexpected value in isimip_simround, filename will be unknown');
end


end % isimip_get_flood_hazard_filename