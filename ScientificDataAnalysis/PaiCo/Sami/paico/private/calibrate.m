function result = calibrate(signal, data, addnoise)

% Shift and scale the signal to the instrumental data
instru = nanmean(data.instrumental.data,1);
origstd = std(signal);

part = signal(ismember(data.target.times, data.instrumental.times));

instru = instru(ismember(data.instrumental.times, data.target.times));
if (nargin == 3 && addnoise)
    if (isfield(data.instrumental,'noisestd') && ~isempty(data.instrumental.noisestd))
        instru = instru + randn(size(instru))*data.instrumental.noisestd;
    end
end

% Method-of-moments aka Mean-variance matching
partmean = mean(part);
instrumean = mean(instru);
mul = std(instru)/std(part);
result.signal = (signal-partmean)*mul + instrumean;

% Ordinary Least Squares with decadal smoothing
% T = zeros(numel(part));
% span = 1;
% T(sub2ind(size(T),floor((0:(numel(part)-1))/span)+1,1:numel(part))) = 1;
% P = [ones(numel(part),1) part(:)];
% 
% alpha = (P'*(T'*T)*P)\(P'*(T'*T)*instru');
% result.signal = alpha(1) + alpha(2)*signal;

result.times = data.target.times;
result.noisestd = std(result.signal)/origstd;
result.uncalibSignal = signal;
