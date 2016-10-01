function [tc_track,tc_track_raw]=isimip_tc_track_load(track_filename,check_plot)
% climada template
% MODULE:
%   module name
% NAME:
%   isimip_tc_track_load
% PURPOSE:
%   load .mat file with (Kerry Emanuel) tracks and convert to climada
%   tc_track strcuture for use in all climada TC functions
%
%   data reference: tobias.geiger@pik-potsdam.de
%
%   next call: climada_tc_hazard_set (or climada_tc_hazard_set_for)
%   possibly also: climada_tc_random_walk...
% CALLING SEQUENCE:
%   tc_track=isimip_tc_track_load(track_filename,check_plot)
% EXAMPLE:
%   tc_track=isimip_tc_track_load('temp_mpi20thcal',0)
%   tc_track=isimip_tc_track_load('temp_mpircp85cal',0)
% INPUTS:
%   track_filename: filename of the .mat file with the (Kerry Emanuel)
%       tracks, default folder is ../climada_data/isimip
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plot: whether show a check plot (=1), or not (=0, default)
%       Note that plotting might often take longer than the full
%       conversion...
% OUTPUTS:
%   tc_track: a tc_track structure, with fields
%   tc_track_raw: the raw data, inspect the structure
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160928, initial
% David N. Bresch, david.bresch@gmail.com, 20160929, renamed to isimip_tc_track_load (from isimip_load_tc_tracks)
%-

tc_track=[];tc_track_raw=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('track_filename','var'),track_filename='';end
if ~exist('check_plot','var'),check_plot=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end 


% template to prompt for track_filename if not given
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
if isempty(fP),fP=isimip_data_dir;end
if isempty(fE),fE='.mat';end
track_filename=[fP filesep fN fE];
if ~exist(track_filename,'file')
    fprintf('Error: %s not found\n',track_filename);
end

% the filanem to save the climada-style tc track to
save_filename=[fP filesep fN '_tc_track' fE];

fprintf('loading from %s ...',track_filename);
tc_track_raw=load(track_filename,'latstore','longstore','vstore','rmstore','pstore','yearstore','monthstore','daystore','hourstore');
tc_track_raw.filename=track_filename;
tc_track_raw.rmstore  =(tc_track_raw.rmstore/1.852); % km to nautical miles
fprintf('done\n');

% figure number of tracks and 'years' represented
n_tracks=size(tc_track_raw.latstore,1);
n_years=n_tracks/600;
fprintf('loaded %i tracks (each %i nodes), representing %i years\n',...
    n_tracks,size(tc_track_raw.latstore,2),n_years);

% 	latsall=(latsall*10).astype('int')
% 	if (hemisphere=='N'):
% 		latsall=latsall[(np.where(latsall > 0))]
% 	else:
% 		latsall=latsall[(np.where(latsall < 0))]
%
% 	if (hemisphere=='N') & ((latsall > 0).any()):
% 		pass
% 	elif (hemisphere=='S') & ((latsall < 0).any()):
% 		#here=True
% 		pass

% 	if ((lonsall > 0).any()) & ((lonsall < 0).any()):
% 		lonspos=lonsall[lonsall >= 0]
% 		lonsneg=lonsall[lonsall < 0]
% 		#here=True
% 		if (((lonspos <= lonmax).any()) & ((lonspos >= lonmin).any())) & (((lonsneg <= lonmax).any()) & ((lonsneg >= lonmin).any())):
% 			pass
% 		else:
% 			stop=True
% 	else:
% 		#print ((lonsall >= lonmin).any())
% 		#print np.min(lonsall) , np.max(lonsall)
% 		if ((lonsall >= lonmin).any()) & ((lonsall <= lonmax).any()):
% 			here=True
% 			pass
% 		else:
% 			stop=True

% for-loop progress to stdout
t0       = clock;
mod_step = 100; % first time estimate after 10 events, then every 100
format_str='%s';

next_track=1;
for track_i=1:n_tracks
    pos=find(abs(tc_track_raw.latstore(track_i,:))>0 & abs(tc_track_raw.longstore(track_i,:))>0);
    if ~isempty(pos)
        
        if min(tc_track_raw.latstore(track_i,pos))>0 % NH only (for TEST)
            
            % store into tc_track
            tc_track(next_track).lon             =tc_track_raw.longstore(track_i,pos);
            tc_track(next_track).lat             =tc_track_raw.latstore(track_i,pos);
            tc_track(next_track).CentralPressure =tc_track_raw.pstore(track_i,pos);
            tc_track(next_track).MaxSustainedWind=tc_track_raw.vstore(track_i,pos);
            tc_track(next_track).RadiusMaxWind   =tc_track_raw.rmstore(track_i,pos);
            
            % convert lon to climada conventino -180..180 instead of 0..360
            lon_pos=find(tc_track(next_track).lon>260);
            if ~isempty(lon_pos),tc_track(next_track).lon(lon_pos)=tc_track(next_track).lon(lon_pos)-360;end
            
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
            tc_track(next_track).datenum=datenum(yyyy,mm,dd,hh,0,0);
          
            %             year_i=climada_global.present_reference_year-n_years+floor(track_i/600);
            %             tc_track(next_track).datenum=datenum(year_i,1,1)+(0:length(tc_track(next_track).lon)-1)*node_timestep/24;
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
    
    % the progress management
    if mod(track_i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/track_i;
        tracks_remaining  = n_tracks-track_i;
        t_projected_sec   = t_elapsed_event*tracks_remaining;
        if t_projected_sec<60
            msgstr = sprintf('converting ... est. %3.0f sec left (%i/%i tracks)',t_projected_sec,   track_i,n_tracks);
        else
            msgstr = sprintf('converting ... est. %3.1f min left (%i/%i tracks)',t_projected_sec/60,track_i,n_tracks);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
    
end % track_i

fprintf(format_str,''); % move carriage to begin of line

if check_plot
    fprintf('plotting %i tracks ...',length(tc_track))
    for track_i=1:length(tc_track)
        plot(tc_track(track_i).lon,tc_track(track_i).lat); hold on
    end
    climada_plot_world_borders
    fprintf(' done\n')
end % check_plot

fprintf('saving tc_track as %s\n',save_filename)
save(save_filename,'tc_track');

end % isimip_tc_track_load