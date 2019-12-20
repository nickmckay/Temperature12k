function result = mom(data, options, tracker)
%MOM Multiproxy climate reconstruction with method of moments
%   result = MOM(data, options, tracker) calculates the climate signal
%   using the method of moments.
%
%   If options is omitted, defaults will be used. If tracker is omitted, no
%   output is generated.
%   
%   If data is set to 'info', result returns default options and
%   information about their meaning.
%
%   For information about data and result structures, type:
%   help rtstruct 

%   COPYRIGHT 
%       2012, Sami Hanhijärvi <sami.hanhijarvi@iki.fi>
%
%   LICENSE
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
%
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
%
%       You should have received a copy of the GNU General Public License
%       along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin == 0
    error('RECON:MOM','At least one input parameter has to be defined.');    
end

if (ischar(data) && strcmp(data,'info'))
    % User requested

    result.info.method.name = 'Method of moments (CPS)';
    result.info.method.description = ...
        ['Composite-plus-scale using method of moments.'];
    result.info.method.authors = 'Various.';
%     result.info.method.reference = ['Lee  TCK,  Zwiers  FW,  Tsao  M  (2007)' ...
%         ' Evaluation  of proxy-based millennial reconstruction methods. ' ...
%         'Climate Dynamics  31(2-3):263–281,  DOI  10.1007/s00382-007-0351-9'];
%     result.info.method.doi = '10.1198/jasa.2010.ap09379';    
    result.options = [];
    return;
end

% Check data
rtstruct(data);

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

%% Combine data to a single matrix
nrdata = 0;
for i = 1:numel(data.proxy)
    nrdata = nrdata + size(data.proxy{i}.data,1);
end
nrtime = numel(data.target.times);
X = nan(nrtime, nrdata);
I = nan(nrtime, 1);
[timemask loc] = ismember(data.instrumental.times, data.target.times);
I(loc(loc > 0),1) = data.instrumental.data(timemask);

pos = 0;
for i = 1:numel(data.proxy)
    if (isfield(data.proxy{i},'lower') && ~isequal(data.proxy{i}.lower, data.proxy{i}.times))
        % This proxy is not annually resolved. We need to interpolate
        for j = 1:size(data.proxy{i}.data,1)
            pos = pos + 1;
            X(:,pos) = interp1(data.proxy{i}.times, data.proxy{i}.data(j,:),data.target.times);
        end
    else
        [mask,loc] = ismember(data.target.times, data.proxy{i}.times);        
        X(mask,pos+(1:size(data.proxy{i}.data,1))) = data.proxy{i}.data(:,loc(loc>0))';
        pos = pos + size(data.proxy{i}.data,1);
    end
end

%% Remove any year that does not have data.
mask = ~all(isnan([I X]),2);
X = X(mask,:);
I = I(mask);
result.times = data.target.times(mask);
% Remove proxies that do not contain data.
mask = ~all(isnan(X),1);
X = X(:,mask);

%% Normalize proxies over common time
mask = all(~isnan(X),2);
X = (X-repmat(mean(X(mask,:),1),size(X,1),1))./repmat(std(X(mask,:),[],1),size(X,1),1);




%% Calculate anomaly and shift & scale it to instrumental data
commonTime = all(~isnan([I X]),2);
if (nnz(commonTime) < 3)
    error('RECON:MOM','Common time period for proxies and instrumental is < 3.');
end
anomaly = nanmean(X,2)';
anomaly = (anomaly-mean(anomaly(commonTime)))/std(anomaly(commonTime));

result.signal = anomaly*std(I(commonTime)) + mean(I(commonTime));