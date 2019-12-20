function argout = varcache(operation, id, data)
% VARCACHE - Single write / multiple read cache
%  ID = VARCACHE(0,file)    Creates a new cache and returns index ID to it. If file is defined, creates or opens the store in that filename.
%       VARCACHE(1,ID,VAR)  Stores VAR to the cache file.
%       VARCACHE(2,ID)      Rewinds the cache ID to the first position.
% VAR = VARCACHE(3,ID)      Retrieves the next variable from cache ID.
%       VARCACHE(4,ID)      Closes the cache ID and destroys it if it was temporary.
%	  	VARCACHE(5)	     Closes all caches and destroys all temporary files associated to it. Does not delete files that were given names when initializing.
%
% Copyright 2011 Sami Hanhijärvi

persistent cache;

types = {'double' 'single' 'int32' 'uint32'};

if (nargin < 1)
    error('At operation number needs to be set.');
end

switch (operation)
    case 0 		
        if (exist('id','var'))
            filename = id;
        else
            filename = tempname;
        end
        if (exist('data','var'))
            opt = data;
        else
            opt = 'w+';
        end
        id = [];
        if (numel(cache) > 0)
            id = find(isnan(cache),1,'first');
        end
        if (isempty(id))
            id = numel(cache)+1;
        end
        cache(id) = fopen(filename, opt);
        argout = id;
		
    case 1
        if (id < 0 || numel(cache) < id || isnan(cache(id)))
            error(['Cache id ' num2str(id) ' not open.']);
        end
        if (~exist('data','var'))
            error('Data must be given');
        end
        typeid = find(strcmp(types, class(data)));
        if (isempty(typeid))
            error(['Data type ' class(data) ' not supported.']);
        end
        
        fwrite(cache(id), ndims(data), 'uint');
        fwrite(cache(id), size(data), 'uint');
        fwrite(cache(id), typeid, 'uint');
        fwrite(cache(id), data, class(data));
    case 2
        if (id < 0 || numel(cache) < id || isnan(cache(id)))
            error(['Cache id ' num2str(id) ' not open.']);
        end
        frewind(cache(id));
    case 3
        if (id < 0 || numel(cache) < id || isnan(cache(id)))
            error(['Cache id ' num2str(id) ' not open.']);
        end
        datadim = fread(cache(id), 1, 'uint');
        datasize = fread(cache(id), datadim, 'uint');
        datasize = (datasize(:))';
        typeid = fread(cache(id), 1, 'uint');
        
        argout = reshape(fread(cache(id),prod(datasize), types{typeid}),datasize);                
    case 4
        for i = 1:numel(id)
            if (id(i) < 0 || numel(cache) < id(i) || isnan(cache(id(i))))
                continue;
            end   
            fclose(cache(id(i)));
            cache(id(i)) = nan;
        end
    case 5
        for i = 1:numel(id)
            if (id(i) < 0 || numel(cache) < id(i) || isnan(cache(id(i))))
                continue;
            end   
            filename = fopen(cache(id(i)));
            fclose(cache(id(i)));
            delete(filename);
            cache(id(i)) = nan;
        end
    otherwise
        error('No such operation');
end


