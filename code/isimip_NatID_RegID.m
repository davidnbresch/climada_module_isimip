function NatID_RegID=isimip_NatID_RegID(ISO3_code,isimip_NatID)
% climada template
% MODULE:
%   isimip
% NAME:
%   isimip_NatID_RegID
% PURPOSE:
%   return the isimip standard list of ISO3 country codes and the isimip
%   internal country number (NatID) as well as the isimip regions for
%   different purposes. 
%
%   Programmer's note: replaces isimip_ISO3_list
%
% CALLING SEQUENCE:
%   ISO3_list=isimip_NatID_RegID(ISO3_code,isimip_NatID)
% EXAMPLE:
%   ISO3_list=isimip_NatID_RegID; % return whole list
%   NatID_RegID_DEU=isimip_NatID_RegID('DEU'); % entries for DEU
%   NatID_RegID_USA_DEU=isimip_NatID_RegID('',[217 52]); % return entries for NatID 52
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   ISO3_code: ISO3 country code(s)
%   isimip_NatID: isimip internal country number(s), so-called NatID
%       if you provide a vector, such as [217 52], get a strcuture with all
%       cuntries that matched (in this example USA and DEU).
% OUTPUTS:
%   NatID_RegID: a strcuture with fields
%       ISO3{i}: the ISO3 country code for counry i
%       NatID(i): the isimip NatID code for counry i
%       TCRegID(i): the isimip region ID for TC
%       TCRegName{i}: the isimip name for the TC region
%   future: further regions, such as FLRegID, FLRegName (for flood), just
%       add the columns in the .xlsx file, the code automatically reads all
%       columns.
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170224, initial, replaces isimip_ISO3_list
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end

% PARAMETERS
%
NatID_RegID_filename=[isimip_data_dir filesep 'NatID_RegID.xlsx']

[ok,NatID_RegID_filename_mat]=climada_check_matfile(NatID_RegID_filename);
if ok
    load(NatID_RegID_filename_mat)
else
    NatID_RegID=climada_xlsread(0,NatID_RegID_filename,'NatD',1);
    % transpose all (since all arrays in climada are 1xn if possible, since
    % we make use of climada_subarray below
    field_names=fieldnames(NatID_RegID);
    for field_i=1:length(field_names)
        if ~ischar(NatID_RegID.(field_names{field_i}))
            NatID_RegID.(field_names{field_i})=NatID_RegID.(field_names{field_i})';
        end
    end % field_i
    save(NatID_RegID_filename_mat,'NatID_RegID');
end


if ~exist('ISO3_code','var'),ISO3_code = '';end
if ~exist('isimip_NatID','var'),isimip_NatID = [];end

if ~isempty(ISO3_code)
    iso3_pos=strmatch(ISO3_code,NatID_RegID.ISO3);
    if ~isempty(iso3_pos)
        NatID_RegID=climada_subarray(NatID_RegID,iso3_pos);
    else
        NatID_RegID=[];
    end % ~isempty(iso3_pos)
    return
end

if ~isempty(isimip_NatID)
    NatID_pos=ismember(NatID_RegID.NatID,isimip_NatID);
    if sum(NatID_pos)>0
        NatID_RegID=climada_subarray(NatID_RegID,NatID_pos);
    else
        NatID_RegID=[];
    end % sum(NatID_pos)>0
    return
end

end % isimip_NatID_RegID