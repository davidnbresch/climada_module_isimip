function entity = fix_coords_isimip(entity,res)
% climada isimip entity fix coordinate shift
% MODULE:
%   isimip
% NAME:
%   fix_coords_isimip
% PURPOSE:
%   Fix the coordinate error in ISIMIP flood and asset data
%
%   Corrects lon_orig and lat_orig to make them match a given resolution
%   (currently '0360as' or '0150as'). This does not do interpolation, just
%   fixes the coordinates that in some cases have very small shifts in
%   different isimip files (e.g. flood depth versus assets).
%
%
%   next call: isimip...
% CALLING SEQUENCE:
%   entity=fix_coords_isimip(entity_orig,res)
% EXAMPLE:
%   clear params;params.val_filename='0150as';
%   [entity,params]=isimip_gdp_entity('DEU',params);
%   entity=fix_coords_isimip(entity, '0150as');
% INPUTS:
%   entity: the original entity whosecoordinates require correction.
% OPTIONAL INPUT PARAMETERS:
%   res: the resolution of the data, currently either '0150as' (default) or
%   '0360as'.
% OUTPUTS:
%   entity: the entity with corrected coordinates (entity.assets.lon and
%   entity.assets.lat)
%
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180301, initial
% David N. Bresch, dbresch@ethz.ch, 20180306, Warning: sizes of lat_true and lat_unique do not match
% David N. Bresch, dbresch@ethz.ch, 20180314, ERROR and WARNING uppercase
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180322, fix for cases
%   where the country is not continuous in lon or lat (e.g. Angola, where a
%   gap between two parts of the country caused an error).
% David N. Bresch, dbresch@ethz.ch, 20180323, \nERROR, WARNING without \n in front
% David N. Bresch, dbresch@ethz.ch, 20180323, length(lon_unique)>1
%-

if ~exist('res','var'); res   = '0150as'; end

lon_orig = entity.assets.lon;
lat_orig = entity.assets.lat;
lon_unique = unique(lon_orig);
lat_unique = unique(lat_orig);

% determine dlon and dlat
if strcmp(res, '0150as')
    dlon = 360/8640;
    dlat   = 180/4320;
elseif strcmp(res, '0360as')
    dlon = 0.1;
    dlat = 0.1;
else
    fprintf('\nERROR: unexpected value in res, aborted\n');
    return
end

% check that dlon,dlat correspond roughly to the entity
if length(lon_unique)>1
    lon_ratio = dlon/(min(diff(lon_unique)));
    if lon_ratio<0.99 || lon_ratio>1.01
        fprintf('WARNING: available and desired resolution might differ (lon), please check\n');
    end
    lat_ratio = dlat/(min(diff(lat_unique)));
    if lat_ratio<0.99 || lat_ratio>1.01
        fprintf('WARNING: available and desired resolution might differ (lat), please check\n');
    end
end
clear lon_ratio lat_ratio;

% get true lon/lat subset
lon_true = (-180+dlon/2):dlon:(180-dlon/2);
lon_true = lon_true(lon_true > min(lon_orig)-dlon/2 & lon_true < max(lon_orig)+dlon/2);
if length(lon_true) ~= length(lon_unique)
    fprintf('WARNING: sizes of lon_true and lon_unique do not match, assuming the country is not continuous in longitude - if not please check\n');
%     return
end
lat_true = (-90+dlat/2):dlat:(90-dlat/2);
lat_true = lat_true(lat_true > min(lat_orig)-dlat/4 & lat_true < max(lat_orig)+dlat/4);
if length(lat_true) ~= length(lat_unique)
    if length(lat_true) == length(lat_unique)-1
        lat_true = (-90+dlat):dlat:(90-dlat);
        lat_true = lat_true(lat_true > min(lat_orig)-dlat/4 & lat_true < max(lat_orig)+dlat/4);
        fprintf('WARNING: sizes of lat_true and lat_unique do not match, fixed (check)\n');
    else
        fprintf('WARNING: sizes of lat_true and lat_unique do not match, assuming the country is not continuous in latitude - if not please check\n');
%         fprintf('\nERROR: sizes of lat_true and lat_unique do not match, aborted\n');
%         length(lat_true)
%         length(lat_unique)
%         lat_true(1),lat_true(end)
%         lat_unique(1),lat_unique(end)
%         return
    end
end

% replace values by their true values
lon = lon_orig;
for i=1:length(lon_unique)
    [~,i_lon_true] = min(abs(lon_true - lon_unique(i)));
    lon(lon == lon_unique(i)) = lon_true(i_lon_true);
end
% for i=1:length(lon_true)
%     lon(lon == lon_unique(i)) = lon_true(i);
% end
entity.assets.lon = lon;
lat = lat_orig;
for i=1:length(lat_unique)
    [temp i_lat_true] = min(abs(lat_true - lat_unique(i)));
    lat(lat == lat_unique(i)) = lat_true(i_lat_true);
end
entity.assets.lat = lat;

end