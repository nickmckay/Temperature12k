function result = regem(data, options, tracker)
%REGEM Multiproxy climate reconstruction with Regularized 
%   Expectation-Maximization
%   result = regem(data, options, tracker) calculates the climate signal
%   that best matches the pairwise comparisons of all proxy data. 
%
%   If options is omitted, defaults will be used. If tracker is omitted, no
%   output is generated.
%   
%   If data is set to 'info', result returns default options and
%   information about their meaning.
%
%   For information about data and result structures, type:
%   help rtstruct 
% 
%   Method is based on articles:
%   Schneider T (2001) Analysis of incomplete climate data: Estimation
%   of mean values and covariance matrices and imputation of missing values. 
%   Journal of Climate 14:853–871.
%
%   Mann ME, Rutherford S, Wahl E, Ammann C (2007) Robustness of proxy-based 
%   climate ?eld reconstruction methods. Journal of Geophysical Research 
%   112(D12), DOI 10.1029/2006JD008272
%
%   Modified by Sami Hanhijärvi (2011)

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
    error('RECON:REGEM','At least one input parameter has to be defined.');    
end

if (ischar(data) && strcmp(data,'info'))
    % User requested
    [options info] = defaults([]);
    result.options = options;
    result.info = info;
    return;
end

% Check data
rtstruct(data);

% Check options and fill in values that are not defined
if (~exist('options','var'))
    options = [];
end
options = defaults(options);  

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

% Current implementation is optimized for TTLS and no other regression
% method is supported although they appear in the code. 
options.regress = 'ttls';


%% Combine data to a single matrix
nrdata = size(data.instrumental.data,1);
for i = 1:numel(data.proxy)
    nrdata = nrdata + size(data.proxy{i}.data,1);
end

nrtime = numel(data.target.times);
X = nan(nrtime, nrdata);
[mask,loc] = ismember(data.target.times, data.instrumental.times);
X(mask,1:size(data.instrumental.data,1)) = data.instrumental.data(:,loc(mask))';
data.instrumental.timeIndices = loc(mask);
pos = size(data.instrumental.data,1);

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

% Remove any year that does not have data.
mask = ~all(isnan(X),2);
X = X(mask,:);
result.times = data.target.times(mask);
% Remove everything that does not contain data or is of equal value
mask = (~all(isnan(X),1));
for i = 1:size(X,2)
    if (~mask(i))
        continue;
    end
    t = X(~isnan(X(:,i)),i);
    if (all(t(1) == t))
        % Remove proxies that contain only one value.
        mask(i) = false;
    end
end
if (isfield(data.instrumental,'locations'))
    result.locations = data.instrumental.locations(mask(1:size(data.instrumental.locations,1)),:);
end
data.instrumental.includeLocations = find(mask(1:size(data.instrumental.data,1)));
X = X(:,mask);

if (strcmpi(options.regress,'ttls'))
    trunc = options.regpar;
end

%% Preprocess

[n, p]       = size(X);
% number of degrees of freedom for estimation of covariance matrix
dofC         = n - 1;            % use degrees of freedom correction      
% if X is a vector, make sure it is a column vector (a single variable)

% get indices of missing values and initialize matrix of imputed values
missing = isnan(X);
nrmissing = nnz(missing);
if nrmissing == 0
    warning('No missing value flags found.')
    return                                      % no missing values
end
[~,kmis]  = ind2sub([n, p], find(missing));


Xmis = nan(n,p);
Xerr = inf(n,p);

  % initial estimates of missing values
if ~isempty(options.Xmis0)
    % substitute given guesses for missing values
    X(missing)  = options.Xmis0(missing);
    [X, M] = center(X);
else
    [X, M] = center(X);
    X(missing)  = zeros(nrmissing, 1);   % fill missing entries with zeros
end
      
if ~isempty(options.C0) 
    C          = options.C0;
else
    C          = X'*X / dofC;      % initial estimate of covariance matrix
end
  
if strcmpi(options.regress,'ttls')
    if (isempty(trunc))    
        [~,S,~] = svd(X./repmat(diag(C)',size(X,1),1));
        trunc = eigenselect(diag(S).^2);
    end
end

% Find patterns of missing data
[missingPatterns,~,rowToMisPattern] = unique(missing,'rows');
firstPattern = 1;
if (all(~missingPatterns(1,:)))
    % Full of data. Skip this in pttls
    firstPattern = 2;
end
iteration = 0;
rdXmis  = Inf;

tracker('RegEM');
while (iteration < options.maxIterations && rdXmis > options.stagTolerance)
    iteration = iteration + 1;

    % initialize for this iteration ...
    CovRes = zeros(p,p);       % ... residual covariance matrix
    peff_ave = 0;                % ... average effective number of variables 

    % Scale variables to unit variance
    if ~isempty(options.weight)
        dataStd = sqrt(diag(options.weight));
    else
        dataStd = sqrt(diag(C)); 
    end
    dataStd(abs(dataStd) < eps) = 1;   % do not scale constant variables
    X = X ./ repmat(dataStd',n,1);
    C = C ./ (dataStd*dataStd');
    
    if strcmpi(options.regress, 'ttls')
      % compute eigendecomposition of correlation matrix
        [V, d]   = peigs(C, options.neigs);
        peff_ave = (dofC - trunc) * nrmissing;
        Sglobal = V(:,trunc+1:end);
        Sglobal = Sglobal * diag(d(trunc+1:end))* Sglobal';
    end
    
    for patIndex=firstPattern:size(missingPatterns,1)         % cycle over records
        misMask = missingPatterns(patIndex,:);
        rowPattern = rowToMisPattern == patIndex;        
        switch options.regress
            case 'mridge'
                % one multiple ridge regression per record
                [B, S, h, peff]   = mridge(C(kavlr{j},kavlr{j}), ...
                             C(kmisr{j},kmisr{j}), ...
                             C(kavlr{j},kmisr{j}), n-1, optreg);

                peff_ave = peff_ave + peff*pm;   % add up eff. number of variables
                dofS     = dofC - peff;          % residual degrees of freedom

                % inflation of residual covariance matrix
                S        = inflation * S;

                % bias-corrected estimate of standard error in imputed values
                Xerr(j, kmisr{j}) = dofC/dofS * sqrt(diag(S))';

            case 'iridge'
                % one individual ridge regression per missing value in this record
                [B, S, h, peff]   = iridge(C(kavlr{j},kavlr{j}), ...
                     C(kmisr{j},kmisr{j}), ...
                     C(kavlr{j},kmisr{j}), n-1, optreg);

                peff_ave = peff_ave + sum(peff); % add up eff. number of variables
                dofS     = dofC - peff;          % residual degrees of freedom

                % inflation of residual covariance matrix
                S        = inflation * S;

                % bias-corrected estimate of standard error in imputed values
                Xerr(j, kmisr{j}) = ( dofC * sqrt(diag(S)) ./ dofS)';
            case 'ttls'
                % truncated total least squares with fixed truncation parameter
                [B, S]   = pttls(V, d, misMask, trunc, Sglobal);

                dofS     = dofC - trunc;         % residual degrees of freedom

                % inflation of residual covariance matrix
                S        = options.inflation * S;

                % bias-corrected estimate of standard error in imputed values
                Xerr(rowPattern, misMask) = repmat(dofC/dofS * sqrt(diag(S))',nnz(rowPattern),1);
        end

        % missing value estimates
        Xmis(rowPattern, misMask)   = X(rowPattern, ~misMask) * B;

        % add up contribution from residual covariance matrices
        inplaceadd(CovRes,S*nnz(rowPattern),misMask);
    end % loop over missing patterns
    
    % rescale variables to original scaling 
    D = repmat(dataStd',n,1);
    X  = X .* D;
    Xerr = Xerr .* D;
    Xmis = Xmis .* D;
%     C          = C .* repmat(D', p, 1) .* repmat(D, 1, p);
    CovRes = CovRes.*(dataStd*dataStd');

    % rms change of missing values
    dXmis = norm(Xmis(missing) - X(missing)) / sqrt(nrmissing);
    
    % relative change of missing values
    nXmis_pre  = norm(X(missing) + M(kmis)') / sqrt(nrmissing);    
    if nXmis_pre < eps
        rdXmis   = Inf;
    else
        rdXmis   = dXmis / nXmis_pre;
    end

    % update data matrix X
    X(missing)  = Xmis(missing);
    
    % re-center data and update mean
    [X, Mup]   = center(X);                  % re-center data
    M          = M + Mup;                    % updated mean vector

    % update covariance matrix estimate
    C          = (X'*X + CovRes)/dofC; 
    
%     figure(3); clf; plot(X(:,2:end)); hold all;
%     plot(X(:,1), 'linewidth',2,'color',[0.7 0.7 0.7]);
%     drawnow;
%     fprintf('%d: %5.9f\n',iteration, rdXmis);

    tracker('RegEM', iteration, options.maxIterations);
end                                        % EM iteration

% add mean to centered data matrix
X  = X + repmat(M, n, 1);  

result.field = nan(size(data.instrumental.data,1),size(X,1));
result.field(data.instrumental.includeLocations,:) = X(:,1:numel(data.instrumental.includeLocations))';
% result.field = X(:,1:size(data.instrumental.data,1))';


