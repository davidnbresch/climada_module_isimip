function [tc_track,tc_track_raw]=isimip_tc_track_load(track_filename,hemisphere,maxlon,check_plot,rmstore_factor)
% climada template
% MODULE:
%   module name
% NAME:
%   isimip_tc_track_load
% PURPOSE:
%   load .mat file with (Kerry Emanuel) tracks and convert to climada
%   tc_track strcuture for use in all climada TC functions
%
%   NOTE: saving the tc_track structure to a new .mat file is currently
%   disabled, as this just wastes space on disk.
%
%   data reference: tobias.geiger@pik-potsdam.de
%
%   next call: climada_tc_hazard_set (or climada_tc_hazard_set_for)
%   possibly also: climada_tc_random_walk...
% CALLING SEQUENCE:
%   tc_track=isimip_tc_track_load(track_filename,hemisphere,maxlon,check_plot)
% EXAMPLE:
%   tc_track=isimip_tc_track_load('temp_mpi20thcal' ,'N',180,1)
%   tc_track=isimip_tc_track_load('temp_mpircp85cal','S',180,1)
%   tc_track=isimip_tc_track_load('Trial1_GB_dkgfdl_20thcal','S',180,1)
% INPUTS:
%   track_filename: filename of the .mat file with the (Kerry Emanuel)
%       tracks, default folder is ../climada_data/isimip/tc_tracks
%       For later runs, the .mat files are each in a subfolder within
%       tc_tracks, hence the code does also handle this, i.e. for files like
%       ../climada_data/isimip/tc_tracks/AAA/AAA.mat, just enter
%       track_filename as 'AAA'
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   hemisphere: 'N' or 'S', or 'both' (default, speeds up)
%   maxlon: the maximum longitude, either =180 (default) or =360
%       For North Atlantic, for example, 180 is much more suitable.
%       360 makes mainly sense for Pacific (islands).
%       Note that climda prefers to deal with -180..180 whenever possible.
%   check_plot: whether show a check plot (=1), or not (=0, default)
%       Note that plotting might often take longer than the full
%       conversion...
%       =-1 for no plot and no messages to stdout (silent)
%   rmstore_factor: factor to multiply rmstore from Kerry's file when
%       storing into tc_track.RadiusMaxWind. Default=1.0
%       Some data files from Kerry reuqire this correction (factor 2), see
%       also rmstore_factor_filelist in PARAMETERS section of this code.
%       To check, you can compare tc_track_raw.rmstore (uncorrected) with
%       tc_track.RadiusMaxWind (correction applied).
% OUTPUTS:
%   tc_track: a tc_track structure, with fields
%   tc_track_raw: the raw data, inspect the structure
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160928, initial
% David N. Bresch, david.bresch@gmail.com, 20160929, renamed to isimip_tc_track_load (from isimip_load_tc_tracks)
% David N. Bresch, david.bresch@gmail.com, 20170402, check_plot=-1
% David N. Bresch, david.bresch@gmail.com, 20170704, rmstore_factor added
% David N. Bresch, david.bresch@gmail.com, 20170825, double() for datenum
% David N. Bresch, david.bresch@gmail.com, 20170826, no new .mat file saved
%-

tc_track=[];tc_track_raw=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('track_filename','var'),track_filename='';end
if ~exist('hemisphere','var'),hemisphere='both';end
if ~exist('maxlon','var'),maxlon=180;end
if ~exist('check_plot','var'),check_plot=0;end
if ~exist('rmstore_factor','var'),rmstore_factor=1;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% the number of tracks per year
n_track_per_year=600; % 600
%
% define the default folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end
%
% Kerry track files in this list require a correction. The radius maximum
% wind (rmstore) needs to be multiplied by factor 2 (see rmstore_factor).
% Therefore, the code checks for this files and sets rmstore_factor=2.
rmstore_factor_filelist={
    'temp_ccsm420thcal.mat'
    'temp_ccsm4rcp85_full.mat'
    'temp_gfdl520thcal.mat'
    'temp_gfdl5rcp85cal_full.mat'
    'temp_hadgem20thcal.mat'
    'temp_hadgemrcp85cal_full.mat'
    'temp_miroc20thcal.mat'
    'temp_mirocrcp85cal_full.mat'
    'temp_mpi20thcal.mat'
    'temp_mpircp85cal_full.mat'
    'temp_mri20thcal.mat'
    'temp_mrircp85cal_full.mat'
    };

% prompt for track_filename if not given
if isempty(track_filename) % local GUI
    track_filename=[isimip_data_dir filesep '*.mat'];
    [filename, pathname] = uigetfile(track_filename, 'Open:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        track_filename=fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(track_filename);
if isempty(fP),fP=[isimip_data_dir filesep 'tc_tracks'];end
if isempty(fE),fE='.mat';end
track_filename=[fP filesep fN fE];
if ~exist(track_filename,'file') % try case of tracks in sub-folders
    track_filename=[fP filesep fN filesep fN fE]; % fN filesep fN is correct
end
if ~exist(track_filename,'file') % does really not exist
    fprintf('Error: %s not found\n',track_filename);
end

% the filanem to save the climada-style tc track to
save_filename=[fP filesep fN '_tc_track_' hemisphere fE];

if exist(save_filename,'file')
    fprintf('loading from previously processed %s ...\n',save_filename);
    load(save_filename)
else
    fprintf('loading from %s ...',track_filename);
    tc_track_raw=load(track_filename,'latstore','longstore','vstore','rmstore','pstore','yearstore','monthstore','daystore','hourstore');
    tc_track_raw.filename=track_filename;
    tc_track_raw.rmstore  =(tc_track_raw.rmstore/1.852); % km to nautical miles
    fprintf('done\n');
    
    % check for RadiusMaxWind
    [~,fN,fE]=fileparts(track_filename);
    if sum(strcmp([fN fE],rmstore_factor_filelist))>0
        rmstore_factor=2;
        fprintf('NOTE: rmstore_factor set =%1.1i for %s\n',rmstore_factor,[fN fE]);
    end
    
    % figure number of tracks and 'years' represented
    n_tracks=size(tc_track_raw.latstore,1);
    n_years=n_tracks/n_track_per_year;
    fprintf('loaded %i tracks (each %i nodes), representing %i years\n',...
        n_tracks,size(tc_track_raw.latstore,2),n_years);
    
    if     strcmpi(hemisphere,'N')
        hemisphere_latmin=  0;hemisphere_latmax=90;
    elseif strcmpi(hemisphere,'S')
        hemisphere_latmin=-90;hemisphere_latmax= 0;
    else % both
        hemisphere_latmin=-90;hemisphere_latmax=90;
    end
    
    if check_plot>=0,climada_progress2stdout;end % init, see terminate below
    next_track=1;
    for track_i=1:n_tracks
        pos=find(abs(tc_track_raw.latstore(track_i,:))>0 & abs(tc_track_raw.longstore(track_i,:))>0);
        if ~isempty(pos)
            
            minlat=min(tc_track_raw.latstore(track_i,pos));
            maxlat=max(tc_track_raw.latstore(track_i,pos));
            if minlat>=hemisphere_latmin && maxlat<=hemisphere_latmax
                
                % store into tc_track
                tc_track(next_track).lon             =tc_track_raw.longstore(track_i,pos);
                tc_track(next_track).lat             =tc_track_raw.latstore(track_i,pos);
                tc_track(next_track).CentralPressure =tc_track_raw.pstore(track_i,pos);
                tc_track(next_track).MaxSustainedWind=tc_track_raw.vstore(track_i,pos);
                tc_track(next_track).RadiusMaxWind   =tc_track_raw.rmstore(track_i,pos)*rmstore_factor;
                
                % convert lon to climada convention -180..180 instead of 0..360
                tc_track(next_track).lon=climada_longitude_unify(tc_track(next_track).lon,maxlon,1);
                
                % create date/time
                mm=tc_track_raw.monthstore(track_i,pos);
                dd=tc_track_raw.daystore(track_i,pos);
                hh=tc_track_raw.hourstore(track_i,pos);
                yyyy=repmat(tc_track_raw.yearstore(track_i),size(mm));
                
                % check for month in next year
                diff_pos=find(diff(mm)<0);
                if ~isempty(diff_pos)
                    try
                        yyyy(diff_pos(1)+1:end)=yyyy(1)+1;
                    catch
                        fpritnf('\WARNING: check %i manually for year end change\n\n',next_track);
                    end
                end
                tc_track(next_track).datenum=datenum(double(yyyy),double(mm),double(dd),double(hh),0,0); % double() 20170825
                
                [tc_track(next_track).yyyy,tc_track(next_track).mm,...
                    tc_track(next_track).dd,tc_track(next_track).hh]=...
                    datevec(tc_track(next_track).datenum);
                
                % comlplete fields as needed by climada
                tc_track(next_track).MaxSustainedWindUnit='kn';
                tc_track(next_track).CentralPressureUnit='mb';
                tc_track(next_track).TimeStep=repmat(6,size(tc_track(next_track).MaxSustainedWind));
                
                tc_track(next_track).ID_no=track_i;
                tc_track(next_track).orig_event_flag=1;
                tc_track(next_track).name=num2str(track_i);
                % not set: extratrop: '*************'
                
                if length(tc_track(next_track).lon)>3,next_track=next_track+1;end % point to the next free track
                %if next_track>1000,break; end % % TEST, stop after 100 tracks filled
            end
        end % ~isempty(pos)
        
        if check_plot>=0,climada_progress2stdout(track_i,n_tracks,100,'tracks');end % update
        
    end % track_i
    if check_plot>=0,climada_progress2stdout(0);end % terminate
    
%     fprintf('saving tc_track as %s\n',save_filename)
%     save(save_filename,'tc_track',climada_global.save_file_version) % for HDF5 format (portability)
    
end % save_filename

if check_plot>0
    fprintf('plotting %i tracks ...',length(tc_track))
    for track_i=1:length(tc_track)
        plot(tc_track(track_i).lon,tc_track(track_i).lat); hold on
    end
    climada_plot_world_borders
    fprintf(' done\n')
end % check_plot

end % isimip_tc_track_load