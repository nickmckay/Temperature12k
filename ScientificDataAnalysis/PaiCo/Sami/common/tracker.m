function tracker(name, pos, total)

persistent Tracker_Data;

if isempty(Tracker_Data) || ~isfield(Tracker_Data,'map') || nargin == 0
    Tracker_Data.map = java.util.HashMap;
    Tracker_Data.times = [];
    tic;
    if (nargin == 0)
        return;
    end
end

if (nargin == 1)
	
    if (Tracker_Data.map.containsKey(name))
        index = Tracker_Data.map.get(name);
    else
        Tracker_Data.times = [Tracker_Data.times; 0 0];
        index = size(Tracker_Data.times,1);
        Tracker_Data.map.put(name, index);
    end
    
    Tracker_Data.times(index,:) = [toc, 0];        
    return;
end

index = Tracker_Data.map.get(name);
times = Tracker_Data.times(index,:);

t = toc-times(1);
if ((t-times(2) > 1) || (pos == total))
    disp([name ' (' num2str(pos) '/' num2str(total) ') : elapsed ' num2str(t,'%10.2f') ' remaining ' num2str(t/pos*(total-pos),'%10.2f') ' total ' num2str(t/pos*total,'%10.2f')]);
    drawnow;
    Tracker_Data.times(index,2) = t;
end