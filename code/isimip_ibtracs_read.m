function tc_track=isimip_ibtracs_read(csv_filename,delimiter)
% climada isimip ibtracs read tc
% MODULE:
%   isimip
% NAME:
%   isimip_ibtracs_read
% PURPOSE:
%   read isimip ibtracs tropical cyclone (TC) track data
%   Either a single track or all tracks within a folder
%
%   next call: climada_tc_hazard_set
% CALLING SEQUENCE:
%   tc_track=isimip_ibtracs_read(csv_filename,delimiter)
% EXAMPLE:
%   csv_filename=[climada_global.modules_dir filesep 'isimip' ...
%       filesep 'data' filesep 'isimip' filesep ...
%       'ibtracs_basin-NA_intp-None_1992230N11325.csv']; % Andrew
%   tc_track=isimip_ibtracs_read(csv_filename);
%   climada_tc_track_info(tc_track,1) % check plot
% INPUTS:
%   csv_filename: the filename of an isimip ibtracs tropical cyclone (TC)
%       track data .csv file, OR the folder name containing such single
%       .csv files (processing all .csv files within)
%       > promted for if not given (for single track)
%       ='TEST' to run TEST mode (just one track, Andrew 1992)
% OPTIONAL INPUT PARAMETERS:
%   delimiter: the delimiter, default is climada_global.csv_delimiter
% OUTPUTS:
%   tc_track: a climada TC track structure, see e.g. climada_tc_read_unisys_database
%       plus the fields RadiusMaxWind, EnvironmentalPressure
% David N. Bresch, david.bresch@gmail.com, 20161203, intial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('csv_filename','var'),csv_filename=[];end % OR:
if ~exist('delimiter','var'),delimiter=climada_global.csv_delimiter;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% set TEST track
if strcmp(csv_filename,'TEST')
    csv_filename=[climada_global.modules_dir filesep 'isimip' ...
        filesep 'data' filesep 'isimip' filesep ...
        'ibtracs_basin-NA_intp-None_1992230N11325.csv']; % Andrew
    fprintf('TEST mode, using %s\n',csv_filename)
end

% prompt for csv_filename if not given
if isempty(csv_filename) % local GUI
    csv_filename=[climada_global.data_dir filesep '*.csv'];
    [filename, pathname] = uigetfile(csv_filename, 'Open:');
    if isequal(filename,0) || isequal(pathname,0)
        tc_track=[]; % init output
        return; % cancel
    else
        csv_filename=fullfile(pathname,filename);
    end
end

track_i=1;

if isdir(csv_filename) % figure whether we deal with a folder
    files=dir([csv_filename filesep '*.csv']);
    
    n_files=length(files);
    
    % template for-loop with waitbar or progress to stdout
    t0       = clock;
    mod_step = 2; % first time estimate after 10 events, then every 100 (see below)
    format_str='%s';
    fprintf('processing %i files\n',n_files);
    
    for file_i=1:n_files
        if ~files(file_i).isdir % is a data file
            single_filename=[csv_filename filesep files(file_i).name];
            tc_track(track_i)=isimip_ibtracs_read(single_filename,delimiter);
            track_i=track_i+1; % point to next free track
            
            if mod(file_i,mod_step)==0  % progress management
                mod_step          = 10;
                t_elapsed_files   = etime(clock,t0)/file_i;
                files_remaining  = n_files-file_i;
                t_projected_sec   = t_elapsed_files*files_remaining;
                if t_projected_sec<60
                    msgstr = sprintf('est. %3.0f sec left (%i/%i files)',t_projected_sec,   file_i,n_files);
                else
                    msgstr = sprintf('est. %3.1f min left (%i/%i files)',t_projected_sec/60,file_i,n_files);
                end
                fprintf(format_str,msgstr); % write progress to stdout
                format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end % progress management
            
        end
    end
    fprintf(format_str,''); % move carriage to begin of line
    
else
    tc_track=[]; % init output
    
    res=climada_csvread(csv_filename,delimiter);
    
    if isempty(res),return;end
    
    % store into tc_track structure
    tc_track(track_i).MaxSustainedWindUnit='kn';
    tc_track(track_i).CentralPressureUnit='kn';
    tc_track(track_i).TimeStep=res.tint;
    tc_track(track_i).RadiusMaxWind=res.rmax;
    tc_track(track_i).lon=res.cgps_lon; % or ngps_lon?
    tc_track(track_i).lat=res.cgps_lat; % or ngps_lat?
    tc_track(track_i).MaxSustainedWind=res.vmax;
    tc_track(track_i).CentralPressure=res.pcen;
    tc_track(track_i).EnvironmentalPressure=res.penv;
    tc_track(track_i).EnvironmentalPressure=res.penv;
    tc_track(track_i).orig_event_flag=res.original_data(1);
    tc_track(track_i).name=res.ibtracsID{1};
    
    % fiddle with the data/time
    yyyy=str2double(res.ibtracsID{1}(1:4)); % figure start year
    tc_track(track_i).yyyy=tc_track(track_i).lon*0+yyyy; % init year
    tc_track(track_i).datenum=datenum(yyyy,1,1)+cumsum(tc_track(track_i).TimeStep)/24; % to define a date
    tc_track(track_i).mm=str2num(datestr(tc_track(track_i).datenum,'mm'));
    tc_track(track_i).dd=str2num(datestr(tc_track(track_i).datenum,'dd'));
    tc_track(track_i).hh=str2num(datestr(tc_track(track_i).datenum,'hh'));
    
    % a unique ID
    tc_track(track_i).ID_no=str2double(res.ibtracsID{1}(1:7));
    
end % isdir(csv_filename)

end % isimip_ibtracs_read