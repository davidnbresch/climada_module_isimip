function dist_km=isimip_haversine(lon1,lat1,lon2,lat2)
% climada isimip haversine
% MODULE:
%   isimip
% NAME:
%   climada_template
% PURPOSE:
%   Calculate the great circle distance between two points
%   on the earth (specified in decimal degrees)
%
%   See also climada_nonspheric_distance_m (the present code just to be
%   100% in line with other isimip code)
%
% CALLING SEQUENCE:
%   dist_km=isimip_haversine(lon1,lat1,lon2,lat2)
% EXAMPLE:
%   dist_km=isimip_haversine(lon1,lat1,lon2,lat2)
% INPUTS:
%   lon1,lat1: point 1 (or vector of points)
%   lon2,lat2: point 2 (always one point, if *1 are vectors, the function
%       returns the distances to point 2)
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   dist_km: distance in km
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161204
%-

dist_km=[]; % init output

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('lon1','var'),return;end
if ~exist('lat1','var'),return;end
if ~exist('lon2','var'),return;end
if ~exist('lat2','var'),return;end

if length(lon2)>1,return;end

% convert decimal degrees to radians
lon1=lon1/180*pi;
lon2=lon2/180*pi;
lat1=lat1/180*pi;
lat2=lat2/180*pi;

% haversine formula
dlon = lon2 - lon1;
dlat = lat2 - lat1;
a = sin(dlat/2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon/2).^2;
c = 2 * asin(sqrt(a));

% 6367 km is the radius of the Earth
dist_km = 6367 * c;

end % isimip_haversine