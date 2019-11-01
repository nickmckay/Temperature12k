%function to grid paleoclimate data

function [totalMedian,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(data,latData,lonData,latgrid,longrid)

%inputs
%data: sheet of data with columns of records. Records should already be
%binned to the same timestep and infilled/have NaNs for missing data.
%Records should also already be normalized and averaged by site if
%necessary.
%latData: latitude coordinates of data (column)
%lonData: longitude coordinates of data (column)
%latgrid% e.g., a 5° grid for the high latitudes: %latgrid=57.5:5:87.5;
%longrid: e.g., associated 5° longitude grid:     %longrid=-180:5:180;

%output
%totalMedian = median of all the grids
%gridGroups = index of grids with data
%gidGroupsLat = Latitude of grids with data
%gridGroupsLon = Longitude of grids with data
%gridMean = timeseries of grids with data


%find the nearest coords to lat and lon grids

%find_nearest.m Changed from find_nearest_coords.m because find_nearest_coords was searching each lat and lon independently
[latmatch,lonmatch]=find_nearest([latData lonData],latgrid',longrid');
uniquefinder=latmatch.^2.*lonmatch.^2;
[uLL,uLLiA,uLLiB]=unique(uniquefinder); %unique
toAverage=data; %data to average

%preallocate some things
gridGroups=cell(1,1);
gridMean=nan(size(toAverage,1),length(uLL));
gridGroupsLat=nan(length(uLL),1);
gridGroupsLon=nan(length(uLL),1);

for i = 1:length(uLL)
    gridGroups{i,1}=find(uLLiB==i);
    gridGroupsLat(i,1)=latmatch(uLLiA(i));
    gridGroupsLon(i,1)=lonmatch(uLLiA(i));
    gridMean(:,i)=nanmean(toAverage(:,gridGroups{i,1}),2);
end

%weight the grids by the cosine of latitude
% weightMat=repmat((cosd(gridGroupsLat)./sum(cosd(gridGroupsLat)))',size(toAverage,1),1);
% totalMean=nansum(weightMat.*gridMean,2);

totalMedian = nanmedian(gridMean,2);
end

