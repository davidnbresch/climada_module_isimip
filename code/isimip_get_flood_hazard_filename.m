function hazard_filename=isimip_get_flood_hazard_filename(flood_filename,entity,isimip_simround,years_range,subtract_matsiro,silent_mode)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_get_flood_hazard_filename
% PURPOSE:
%   provide a unique and consistent name for the hazard file (.mat) with
%   flood footprints for given entity, which can be called in scripts and
%   in functions. This avoid using different conventions or repeating code
%   too much.
%
%   next call: isimip_flood_load
% CALLING SEQUENCE:
%   hazard_filename=isimip_get_flood_hazard_filename(flood_filename,entity,check_plot,isimip_data_subdir,years_range)
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
%   isimip_simround: the sub-directory within the isimip folder within
%       climada_data/isimip where the raw isimip data (NetCDF file
%       flood_filename) is located. This does not affect where the hazard is
%       saved. If not specified, the file is assumed to be located directly
%       within the climada_data/isimip folder.
%   years_range: vector of length 2 containing the first and the last year
%       to be loaded from the netcdf file. If empty or [0 0], loads all data.
%   subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%       fraction from the data. Default =0.
%   silent_mode: =1 do not print anything if the file does not exist
% OUTPUTS:
%   hazard_filename: a file name for the .mat hazard file on entity grid.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20171130, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180118, changed a few
% things to be consistent.
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180518, new function
%   argument 'subtract_matsiro'
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180918, adding argument
%   silent_mode
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('flood_filename','var'),         flood_filename=         '';end
if ~exist('entity','var'),                 entity=                 '';end
if ~exist('years_range','var'),            years_range=         [0 0];end
if ~exist('isimip_simround','var'),        isimip_simround=        '';end
if ~exist('subtract_matsiro','var'),       subtract_matsiro=        0;end
if ~exist('silent_mode','var'),            silent_mode=             1;end

% define paths
isimip_data_dir = [climada_global.data_dir filesep 'isimip' filesep isimip_simround];
isimip_hazard_dir = [climada_global.hazards_dir filesep 'isimip' filesep isimip_simround];


% check validity of arguments
if ~isequal(size(years_range), [1 2])
    warning('** error ** years_range should be of size [1 2] *****')
    return
end

% create directory if does not exist
if ~isdir(isimip_hazard_dir)
    mkdir(isimip_hazard_dir); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_hazard_dir);
end

% prompt for flood_filename if not given
if isempty(flood_filename) % local GUI
    flood_filename=[isimip_data_dir filesep '*.nc'];
    [filename, pathname] = uigetfile(flood_filename, 'Select flood file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        flood_filename=fullfile(pathname,filename);
    end
end

% flood_filename: complete path, if missing
[fP,fN,fE]=fileparts(flood_filename);
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
flood_filename=[fP filesep fN fE];
if ~silent_mode
    if ~exist(flood_filename,'file')
        fprintf('Error: %s not found\n',flood_filename);
    end
end


[~,fN]=fileparts(entity.assets.filename);
[~,fN2]=fileparts(flood_filename);
% fN2 should start with 'flddph', in which case replace with 'FLOOD'
fN2=strrep(fN2, 'flddph', 'FLOOD');
fN=strrep(fN,'_entity','');
if subtract_matsiro fN3='_mFRCmatsiro';, else fN3='';end;
if isequal(years_range, [0 0])
    hazard_filename=[isimip_hazard_dir filesep fN2 '_' fN fN3 '_FL.mat'];
else
    hazard_filename=[isimip_hazard_dir filesep fN2 '_' fN '_' num2str(years_range(1)) '-' num2str(years_range(2)) fN3 '_FL.mat'];
end

% hazard_filename: complete path, if missing
[fP,fN,fE]=fileparts(hazard_filename);
if isempty(fP),fP=climada_global.hazards_dir;end
if isempty(fE),fE='.mat';end
hazard_filename=[fP filesep fN fE];

end % isimip_get_flood_hazard_filename