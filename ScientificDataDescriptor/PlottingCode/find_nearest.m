function [lat,lon] =find_nearest(coords,latvec,lonvec)
%updated and simplified from find_nearest_coords.m. FUnction locates the
%nearest set of coordinates in latvec and lonvec to coords. 

%Inputs: 
%latvec and lonvec repersent the lat and lon vectors of a grid for example. 
%Coords are the coordinates for a subset of data (e.g.) proxy record
%locations that you want to know which gridpoint thay fall in or closest
%to. 

for i = 1:length(coords)
    for j = 1:length(latvec)
        distmat(j,1) = lldistkm(coords(i,:),[latvec(j,1) lonvec(j,1)]);
    end


a=find(distmat == min(distmat));

lat(i)=latvec(a);

lon(i)=lonvec(a);
clear distmat
end


lat=lat';
lon=lon';
end

    
            