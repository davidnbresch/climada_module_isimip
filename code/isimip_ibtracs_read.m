function [tc_track,save_file]=isimip_ibtracs_read(csv_filename,delimiter,save_flag)
% climada isimip ibtracs read tc
% MODULE:
%   isimip
% NAME:
%   isimip_ibtracs_read
% PURPOSE:
%   read isimip ibtracs tropical cyclone (TC) track data
%   Either a single track or all tracks within a folder
%
%   The code basically assumes single .csv files to be stored in
%   {climada_global.data_dir}/isimip/ibtracs/{basin}/*.cs
%   where climada_global.data_dir is the local users climada data dir
%   and basin denotes the ocean basin, such as 'EP','NA','NI','SA','SI','SP','WP'
%
%   If a whoe folder is processed, tracks with less than 3 nodes are
%   skipped.
%
%   next call: climada_tc_hazard_set
% CALLING SEQUENCE:
%   tc_track=isimip_ibtracs_read(csv_filename,delimiter)
% EXAMPLE:
%   tc_track=isimip_ibtracs_read('TEST'); % returns test track (Andrew,1992)
%   tc_track=isimip_ibtracs_read('NA'); % all North Atlantic tracks
%   csv_filename=[climada_global.modules_dir filesep 'isimip' ...
%       filesep 'data' filesep 'isimip' filesep ...
%       'ibtracs_basin-NA_intp-None_1992230N11325.csv']; % Andrew
%   tc_track=isimip_ibtracs_read(csv_filename);
%   climada_tc_track_info(tc_track,1) % check plot
%
% INPUTS:
%   csv_filename: the filename of an isimip ibtracs tropical cyclone (TC)
%       track data .csv file, OR the folder name containing such single
%       .csv files (processing all .csv files within)
%       > promted for if not given (for single track)
%       ='TEST' to run TEST mode (just one track, Andrew 1992)
%        Note that this reads the test file in the isimip module's data folder, i.e
%        {climada_global.modules_dir}/isimip/data/ibtracs/ibtracs_basin-NA_intp-None_1992230N11325.csv
%       ='EP', 'NA', 'NI', 'SA', 'SI', 'SP', 'WP': read all tracks of one
%        basin, assuming the .csv files are named  {climada_global.data_dir}/isimip/ibtracs/ibtracs_basin-{basin}_intp-None_*.csv
%        Resulting tc_track saved as ibtracs{basin}.mat in the same folder the ibtracs file were found
%       ='all' to process all tracks globally in {climada_global.data_dir}/isimip/ibtracs. 
%        Resulting tc_track saved as ibtracs.mat in the same folder the ibtracs file were found
%        Usually processing 6462 single track files, 6323 tracks ok, 139 not ok (less than 3 nodes, skipped)
%        1407 tracks more than once, duplicates removed
% OPTIONAL INPUT PARAMETERS:
%   delimiter: the delimiter, default is climada_global.csv_delimiter
%   save_flag: if =1, save as .mat file, named */ibtracs.mat in the folder
%       that has been processed (only if called for a folder). Default=0
% OUTPUTS:
%   tc_track: a climada TC track structure, see e.g. climada_tc_read_unisys_database
%       plus the fields RadiusMaxWind, EnvironmentalPressure
%   save_file: the file (with path) where the tc_track structure has been
%       saved to if save_flag=1, ='' otherwise
% David N. Bresch, david.bresch@gmail.com, 20161203, intial
% David N. Bresch, david.bresch@gmail.com, 20161222, new field isotime used to properly define yyyy,mm,dd and hh
% David N. Bresch, david.bresch@gmail.com, 20161226, allow for basin name, such as 'NA'
% David N. Bresch, david.bresch@gmail.com, 20170121, allow for all tracks globally (currently basin name not operational), duplicates removed
%-

save_file='';

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('csv_filename','var'),csv_filename=[];end % OR:
if ~exist('delimiter','var'),delimiter='';end
if ~exist('save_flag','var'),save_flag=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
if isempty(delimiter),delimiter=climada_global.csv_delimiter;end
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

basin_name=''; % usually no specific basin

if length(csv_filename)==2 % length =2 for 'XX', such as 'NA'
    % it contains a basin name
    fprintf('processing %s basin\n',csv_filename);
    basin_name=csv_filename;
    csv_filename=[climada_global.data_dir filesep 'isimip' filesep ...
        'ibtracs'];
%        'ibtracs' filesep csv_filename];
    save_flag=1;
end

if strcmpi(csv_filename,'all')
    % it contains a basin name
    fprintf('processing all tracks\n');
    csv_filename=[climada_global.data_dir filesep 'isimip' filesep ...
        'ibtracs'];
    save_flag=1;
end

track_i=1;

if isdir(csv_filename) % figure whether we deal with a folder
    
    if ~isempty(basin_name)
        % files such as e.g. ibtracs_basin-WP_intp-None_2014015N10129
        files=dir([csv_filename filesep 'ibtracs_basin-' basin_name '_intp-None_*.csv']);
    else
        files=dir([csv_filename filesep '*.csv']);
    end
    
    n_files=length(files);
    
    if n_files>0
        
        save_file=[csv_filename filesep 'ibtracs' basin_name '.mat'];
        if exist(save_file,'file')
            fprintf('HINT: consider isimip_ibtracs_load, since %s exists already\n',save_file);
        end
        
        track_ok=0;track_not_ok=0; % init
        
        % template for-loop with waitbar or progress to stdout
        t0       = clock;
        mod_step = 2; % first time estimate after 10 events, then every 100 (see below)
        format_str='%s';
        fprintf('processing %i single track files\n',n_files);
        
        for file_i=1:n_files
            if ~files(file_i).isdir % is a data file
                single_filename=[csv_filename filesep files(file_i).name];
                tc_track(track_i)=isimip_ibtracs_read(single_filename,delimiter);
                
                if length(tc_track(track_i).lon)>2
                    track_i=track_i+1; % point to next free track
                    track_ok=track_ok+1;
                else
                    track_not_ok=track_not_ok+1;
                end
                
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
        
        fprintf('%i tracks ok, %i not ok (less than 3 nodes, skipped)\n',track_ok,track_not_ok);
        
        % check for duplicate entries
        ID_no=zeros(1,length(tc_track));
        for track_i=1:length(tc_track)
            ID_no(track_i)=tc_track(track_i).ID_no;
        end
        [ID_no_unique,unique_i] = unique(ID_no);
        if length(ID_no_unique)~=length(ID_no)
            fprintf('%i tracks more than once, duplicates removed\n',length(ID_no)-length(ID_no_unique));
            tc_track=tc_track(unique_i);
        end
        
        if save_flag
            fprintf('savig tc_track in %s\n',save_file);
            save(save_file,'tc_track');
        end % save_flag
        
    else
        fprintf('WARNING: no .csv files in %s\n',csv_filename);
        tc_track=[];
    end % n_files
    
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
    
    if isfield(res,'isotime')
        isotime=res.isotime;
        tc_track(track_i).yyyy=fix(isotime/1e6);
        isotime               =isotime-tc_track(track_i).yyyy*1e6;
        tc_track(track_i).mm  =fix(isotime/1e4);
        isotime               =isotime-tc_track(track_i).mm  *1e4;
        tc_track(track_i).dd  =fix(isotime/1e2);
        tc_track(track_i).hh  =isotime-tc_track(track_i).dd  *1e2;
        tc_track(track_i).datenum=datenum(tc_track(track_i).yyyy,...
            tc_track(track_i).mm,tc_track(track_i).dd,tc_track(track_i).hh,0,0); % convert
    else
        % fiddle with the data/time
        yyyy=str2double(res.ibtracsID{1}(1:4)); % figure start year
        tc_track(track_i).yyyy=tc_track(track_i).lon*0+yyyy; % init year
        tc_track(track_i).datenum=datenum(yyyy,1,1)+cumsum(tc_track(track_i).TimeStep)/24; % to define a date
        tc_track(track_i).mm=str2num(datestr(tc_track(track_i).datenum,'mm'));
        tc_track(track_i).dd=str2num(datestr(tc_track(track_i).datenum,'dd'));
        tc_track(track_i).hh=str2num(datestr(tc_track(track_i).datenum,'hh'));
    end
    
    % a unique ID
    tc_track(track_i).ID_no=str2double(res.ibtracsID{1}(1:7));
    
    % deal with missing pressure
    tc_track(track_i).CentralPressure(tc_track(track_i).CentralPressure<=0)=NaN;
    if sum(isnan(tc_track.CentralPressure))
        tc_track.CentralPressure=extra_p(tc_track(track_i).MaxSustainedWind,tc_track(track_i).lat,tc_track(track_i).lon);
    end
    
end % isdir(csv_filename)

end % isimip_ibtracs_read


% follow helper functions (the ones easy to convert from python ;-)
% -----------------------------------------------------------------

function p_from_v=extra_p(vmax,lat,lon)
% determined in hPa
%p_from_v = 1024.688+0.055*lat-0.028*lon-0.815*vmax # peduzzi
p_from_v = 1024.388 + 0.047*lat - 0.029*lon - 0.818*vmax; % ibtracs 1980 -2013 (r2=0.91)
end %  extra_p
