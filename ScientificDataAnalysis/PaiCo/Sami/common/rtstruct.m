function data = rtstruct(data)
%RTSTRUCT Recontoolbox structure utility function
%   data = rtstruct() returns an example dataset that can be used with
%   most reconstruction methods. 
%
%   rtstruct(data) checks if the input data conforms to the structure
%   format. Throws an error if it does not.
%   
%   Structure of data: (optional fields are marked with *)
%   data.proxy                      
%       Cell array of proxy data
%   data.proxy{i}.data              
%       Matrix of proxy record. Columns equal times and rows equal proxy 
%       records. Records with identical temporal structure may be stacked.
%       Can contain NaN's.
%   data.proxy{i}.times             
%       A row vector of mid-times for each sample. Must have equal number 
%       of columns as data. Measured in years AD.
%   *data.proxy{i}.lower            
%       Lower, or older, bound in time for each proxy sample. Measured in 
%       years AD.
%   *data.proxy{i}.upper            
%       Upper, or younger, bound in time for each proxy sample. Measured in 
%       years AD.
%   *data.proxy{i}.intargetres      
%       Boolean stating if the record is in the same resolution as the 
%       target. Some methods handle such records differently.
%   *data.proxy{i}.locations        
%       Locations of proxy records. Must have equal number of rows as data.
%       First column is degrees in latitude, second in longitude.
%   data.instrumental              
%       Struct of instrumental data.
%   data.instrumental.data 
%       Instrumental values. Columns correspond to times and rows to
%       locations. Can contain NaN's.
%   data.instrumental.times        
%       Mid-times for samples in instrumental data. Measured in years AD.
%       Currently only annual instrumental data is considered.
%   *data.instrumental.locations    
%       Locations of instrumental data. Must have equal number of rows as
%       data. First column is degrees in latitude, second in longitude.
%   data.target
%       Struct of reconstruction target.
%   data.target.times
%       Mid-times for target of reconstruction. Measured in years AD.
%       Currently only annual targets are considered.

if (nargin == 1)
    assert(isstruct(data), 'Data must be a struct. Type "help rtstruct" for details.');
    
    % Check proxies
    assert(isfield(data,'proxy') && iscell(data.proxy), 'Data must contain a field called proxy that is a a cell array of proxy data. Type "help rtstruct" for details.');
    for i = 1:numel(data.proxy)
        assert(isstruct(data.proxy{i}),['Proxy ' num2str(i) ' is not a struct. Type "help rtstruct" for details.']);
        assert(isfield(data.proxy{i},'data'), ['Proxy ' num2str(i) ' doesn''t contain data field. Type "help rtstruct" for details.']);
        assert(isfield(data.proxy{i},'times'), ['Proxy ' num2str(i) ' doesn''t contain times field. Type "help rtstruct" for details.']);
        assert(size(data.proxy{i}.data,2) == size(data.proxy{i}.times,2), ['Proxy ' num2str(i) ' doesn''t have equal number of columns for data and times fields. Type "help rtstruct" for details.']);        
        if (isfield(data.proxy{i},'lower'))
            assert(isfield(data.proxy{i},'upper'),['Proxy ' num2str(i) ' must contain upper bounds in time if it contains lower bounds. Type "help rtstruct" for details.']);
            assert(size(data.proxy{i}.lower,2) == size(data.proxy{i}.times,2), ['Proxy ' num2str(i) ' doesn''t have equal number of columns for lower and times fields. Type "help rtstruct" for details.']);
            assert(size(data.proxy{i}.upper,2) == size(data.proxy{i}.times,2), ['Proxy ' num2str(i) ' doesn''t have equal number of columns for upper and times fields. Type "help rtstruct" for details.']);            
        else
            assert(~isfield(data.proxy{i},'upper'), ['Proxy ' num2str(i) ' must contain lower bounds in time if it contains upper bounds. Type "help rtstruct" for details.']);
        end
        if (isfield(data.proxy{i},'locations'))
            assert(size(data.proxy{i}.data,1) == size(data.proxy{i}.locations,1),['Proxy ' num2str(i) ' doesn''t have equal number of rows for data and locations. Type "help rtstruct" for details.']);
            assert(size(data.proxy{i}.locations,2) == 2, ['Proxy ' num2str(i) ' doesn''t have two columns in locations. Type "help rtstruct" for details.']);
        end
        if (isfield(data.proxy{i},'intargetres'))
            assert(islogical(data.proxy{i}.intargetres),['Proxy ' num2str(i) ' has a non logical value for field intargetres. Type "help rtstruct" for details.']);
        end
    end
    
    % Check instrumental
    assert(isfield(data,'instrumental') && isstruct(data.instrumental), 'Data must contain a field called instrumental that is a struct of instrumental data. Type "help rtstruct" for details.');
    assert(isfield(data.instrumental,'data'),['Instrumental field doesn''t contain data. Type "help rtstruct" for details.']);
    assert(isfield(data.instrumental,'times'),['Instrumental field doesn''t contain times. Type "help rtstruct" for details.']);
    assert(size(data.instrumental.data,2) == size(data.instrumental.times,2), ['Instrumental field doesn''t have equal numeber of columns for data and times. Type "help rtstruct" for details.']);
    
    % Check target
    assert(isfield(data,'target') && isstruct(data.target), 'Data must contain a field called target that is a struct. Type "help rtstruct" for details.');
    assert(isfield(data.target,'times'), 'Target field must contain a field that contains a field called times. Type "help rtstruct" for details.'); 
    
elseif (nargout == 1)
    % Construct example dataset
    
    % Initial values for parameters
    phi_P = exp(-5);
    tau = 0.1;

    % Proxy parameters
    snr = 5;
    nrproxy = 20;

    % Construct target grid data.
    [lat long] = meshgrid(-20:5:20,-20:5:20);
    data.instrumental.locations = [lat(:) long(:)];
    data.target.locations = data.instrumental.locations;
    data.target.times = 0:100;

    anglescale = 1/2;

    % Create 5 seeds, 1 in each corner and one in the center
    seeds = zeros(5,numel(data.target.times));
    for i = 1:size(seeds,1)
        seeds(i,:) = randn*sin(anglescale * randn*data.target.times / 20*pi+rand) + randn*cos(anglescale * randn*data.target.times/ 20*pi + rand);
    end
    
    % Periodical functions introduce autocorrelation automatically
    seedloc = [min(lat(:)) min(long(:)); max(lat(:)) min(long(:)); ...
        min(lat(:)) max(long(:)); max(lat(:)) max(long(:)); ...
        mean(lat(:)) mean(long(:))];
    D = zeros(size(data.instrumental.locations,1), numel(data.target.times));
    for i = 1:size(seeds,2)
        D(:,i) = griddata(seedloc(:,1), seedloc(:,2), seeds(:,i), lat(:), long(:));
    end
    target = D;
    data.target.data = target;

    % Construct instrumental data
    D = D + std(D(:))*tau*randn(size(D));
    % Instrumental data is only available in the latter 25% of time
    indices = round(numel(data.target.times)*0.75):size(D,2);
    data.instrumental.times = data.target.times(indices);
    data.instrumental.data = D(:,indices);


    % Construct proxy data
    loc = zeros(nrproxy,2);
    dt = zeros(nrproxy, size(target,2));
    for i = 1:nrproxy
        if (i < 5)
            pos = seedloc(i,:);
        else        
            pos = [rand*(seedloc(4,1)-seedloc(1,1))+seedloc(1,1) ...
                rand*(seedloc(4,2)-seedloc(1,2))+seedloc(1,2)];
        end
        dist = earthDistances(pos, data.instrumental.locations);
        S = exp(-phi_P*dist);
        S = S./sum(S);
        t = S*target;
        e = randn(size(t));
        t = t + var(t)/snr*e/var(e);
        loc(i,:) = pos;


        % Random monotonic transformation
        z = randn(size(t));
        [~,ind] = sort(t);
        t(ind) = sort(z);
        t = t + randn(size(t))*std(t)/snr;
        dt(i,:) = (t-mean(t))./std(t);
        % Store data
        data.proxy{i}.locations = pos;
        data.proxy{i}.data = t;
        data.proxy{i}.times = data.target.times;
    end    
else
    eval('help rtstruct');
    return;
end
