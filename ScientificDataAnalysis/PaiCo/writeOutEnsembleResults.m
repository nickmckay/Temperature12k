ar = struct();

latBins = ({'60to90N', '30to60N', '0to30N', '0to30S', '30to60S', '60to90S'});
for j = 1:6
    ar.(['lb' latBins{j}]).mat =  nan(120,500);
    
end


for i = 1:500
    for j = 1:6
        
       ar.(['lb' latBins{j}]).mat(:,i) =  apcb{i,1}(j).signal;  
    
    end
end

for j = 1:6
    csvwrite([latBins{j} '500.csv'],[(11900:-100:0)' ar.(['lb' latBins{j}]).mat])   
end
