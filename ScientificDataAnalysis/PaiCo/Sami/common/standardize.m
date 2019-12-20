function data = standardize(data, method)

assert(nargin == 2, 'Two arguments required.');

switch lower(method)
    case {'common'}
        % Check conditions
        assert(size(data.instrumental,2) < 2, ...
            ['Common time calibration can only be carried out if the ' ...
             'instrumental record is a single global signal.']);
         
        times = data.instrumental.times;
        for i = 1:numel(data.proxy)
            times = intersect(times, data.proxy{i}.times);
        end

        mask = ismember(data.instrumental.times, times);
        X = data.instrumental.data;
        X = (X-repmat(nanmean(X(:,mask),2),1,size(X,2)))./repmat(nanstd(X(:,mask),[],2),1,size(X,2));
        data.instrumental.data = X;

        for i = 1:numel(data.proxy)
            X = data.proxy{i}.data;
            X = (X-repmat(nanmean(X(:,mask),2),1,size(X,2)))./repmat(nanstd(X(:,mask),[],2),1,size(X,2));
            data.proxy{i}.data = X;
        end        
    case 'individual'        
        for i = 1:numel(data.proxy)
            X = data.proxy{i}.data;
            X = (X-repmat(nanmean(X,2),1,size(X,2)))./repmat(nanstd(X,[],2),1,size(X,2));
            data.proxy{i}.data = X;
        end
    case 'closest'
        
    
        
    otherwise 
        error(['Method "' method '" not recognized.']);
end