function hazard=ncar_tc_wind_import(filename,misdat_value,min_value)
% climada NCAR TC windfield import
% MODULE:
%   isimip
% NAME:
%   ncar_tc_wind_import
% PURPOSE:
%   read NCAR single event tropical cyclone (TC) footprint files.
%
%   Format:
%   wind(m/s),    lat,     lon
%   -99.00,   21.79,  270.90
%   -99.00,   21.79,  270.94
%
%   CAUTION: does NOT YET handle windfields on different grids or different
%   sub-areas of a grid.
%
%   next call: climada_hazard_plot
% CALLING SEQUENCE:
%   hazard=ncar_tc_wind_import(filename)
% EXAMPLE:
%   hazard=ncar_tc_wind_import; % TEST, reads demo windfield
%   p.max_value=80;climada_hazard_plot_nogrid(hazard,1,2,p); % show windfield
% INPUTS:
%   filename: filename of a NCAR TC wind file or a folder, in which case
%       all files within are imported
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   misdat_value: the value for missing data, default =-99.00
%   min_value: minimum value to keep (since very low windspeeds do not
%       result in any impact). Default=10 m/s, values<min_value are set =0
% OUTPUTS:
%   hazard: a climada hazard event set, see manual for all fields, core
%       fields of this struct are:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171020, initial
%-

hazard=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables


% poor man's version to check arguments
if ~exist('filename','var'),filename=[];end
if ~exist('misdat_value','var'),misdat_value=[];end
if ~exist('min_value','var'),min_value=[];end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% set default value for filename (for TEST mode)
if isempty(filename),filename=[module_data_dir filesep 'isimip' filesep 'NCAR_frances_200490418_mod.txt'];end
%
% default missing data value
if isempty(misdat_value),misdat_value=-99.00;end
if isempty(min_value),min_value=10;end % m/s
    
if ~exist(filename,'file')
    fprintf('ERROR: file not found, aborted: %s\n',filename);
    return
elseif isdir(filename)
    files=dir(filename);
    filenames={}; % init
    for file_i=1:length(files)
        if ~files(file_i).isdir && ~strcmpi(files(file_i).name(1),'.')
            filenames{end+1}=[filename filesep files(file_i).name];
        end
    end %file_i
else
    filenames{1}=filename;
end

for file_i=1:length(filenames)
    data=importdata(filenames{file_i},',',1);
    hazard.lon = data.data(:,3)';
    hazard.lat = data.data(:,2)';
    data.data(data.data(:,1)==misdat_value,1)=0; % replace misdat by zero
    data.data(data.data(:,1)<min_value,1)=0; % replace misdat by zero
    hazard.intensity=sparse(data.data(:,1)');
end % file_i

% make sure longitudes are canonical
hazard.lon(hazard.lon>180)  = hazard.lon(hazard.lon>180)-360;
hazard.lon(hazard.lon<-180) = hazard.lon(hazard.lon<-180)+360;

% complete the hazard structure
hazard.centroid_ID=1:length(hazard.lon);
hazard.peril_ID='TC';
hazard.date=datestr(now);
hazard.comment=sprintf('TC hazard event set, generated %s by %s',hazard.date,mfilename);
hazard.orig_years=NaN;
hazard.event_count=size(hazard.intensity,1);
hazard.orig_event_count=hazard.event_count;
hazard.event_ID=1:hazard.event_count;
hazard.orig_event_flag=hazard.event_ID*0+1;
hazard.frequency=hazard.orig_event_flag;
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
hazard.windfield_comment=filename;
hazard.filename=filename;
hazard.reference_year=climada_global.present_reference_year;
hazard.fraction=spones(hazard.intensity); % fraction 100%

end % ncar_tc_wind_import