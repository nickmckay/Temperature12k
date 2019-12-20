function monotonicity

darkblue = [42 75 92]/255;
lightblue = [40 175 245]/255;
cyan = [12 194 158]/255;
brown = [206 82 64]/255;
red = [245 40 57]/255;

% darkblue = [29 76 76]/255;
% lightblue = [50 65 255]/255;
% cyan = [73 178 40]/255;
% brown = [75 33 76]/255;
% red = [255 54 35]/255;


x = linspace(0,5,1e2);
y = zscore(sinc(x));
Y = repmat(y,10,1);
Y = Y + randn(size(Y))*0.2;

X = zeros(size(Y));
for i = 1:size(X,1)
    f = zscore(exp(zscore(sort(gamrnd(rand*10, rand*10,1,numel(x))))));
    [~,ind] = sort(Y(i,:));
    X(i,ind) = zscore(f);
end

data.proxy{1}.data = X;
data.proxy{1}.lower = 1:numel(x);
data.proxy{1}.upper = (1:numel(x))+1;
data.instrumental.times = 1:numel(x);
data.instrumental.data = randn(1,numel(x));
data.target.times = 1:numel(x);
options.autocorrelation = 0;

result = paico('run', data, options, @tracker);

figure(1); clf; hold on;
plot(x,y,'-','color',red,'linewidth',1);
plot(x,mean(X,1),':','color',lightblue,'linewidth',1);
plot(x,zscore(result.signal), '-','color',darkblue,'linewidth',1);
legend({'sinc' 'CPS' 'PaiCo'});
squarepage([8,4]);
print('-depsc2','transfer.eps');

figure(2); clf;
plot(X');

%     yt = y;
%     yt(ind) = f;
%     figure(1); clf; 
%     subplot(1,2,1); hold on;   
%     plot(x,val,'r-');
%     plot(x,f,'k-');
%     subplot(1,2,2); hold on;
%     plot(x,y,'r-');
%     plot(x,yt,'k-');
%     pause;
