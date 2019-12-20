function distances = earthDistances(rowloc, colloc)
%% Calculate distances matrix between geographical locations
%Input is N by 2, each row a (lat, long) pair, -90<lat<90; -180<long<180.
%Output is a N by N matrix of great circle distances in KM (approximating the
%earth as a sphere), where the (i,j) entry is the distance between the ith
%and jth rows of the input vector. So the diagonal is zero. 
%This makes use of the so-called haversine formulation (see wikipedia),
%which is also used in the m_lldist.m code of the m_map package. (m_lldist
%gives identical results, but doesn't seem well set up for the formulation
%of the matrix we want here.)
%radius of the earth si taken as 6378.137

R = 6378.137; %radius of the earth in km

if (~exist('colloc','var'))
    colloc = rowloc;
end

rowloc = rowloc/180*pi;
colloc = colloc/180*pi;
distances = zeros(size(rowloc,1),size(colloc,1));
for i = 1:size(colloc,1)
    distances(:,i) = R*2*asin(sqrt(sin((rowloc(:,1)-colloc(i,1))/2).^2 ...
        + cos(rowloc(:,1)).*cos(colloc(i,1)).*sin((rowloc(:,2)-colloc(i,2))/2).^2));
end