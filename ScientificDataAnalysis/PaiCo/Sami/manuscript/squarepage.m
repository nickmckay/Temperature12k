function squarepage(dims)

if (nargin < 1)
    dims = [8.2 8.2];
end

smallfontsize = 7;
fontsize = 8;
axs = findobj(gcf,'type','axes');
for h = axs
    set(h,'FontSize',fontsize);
    set(findobj(h,'Type','text'),'FontSize',smallfontsize);
    tls = get(h,'title');
    for t = 1:length(tls)
        if (iscell(tls))
            set(tls{t},'FontSize',fontsize);
        else 
            set(tls,'FontSize',fontsize);
        end
    end
    tls = get(h,'xlabel');
    for t = 1:length(tls)
        if (iscell(tls))
            set(tls{t},'FontSize',fontsize);
        else 
            set(tls,'FontSize',fontsize);
        end
    end
    tls = get(h,'ylabel');
    for t = 1:length(tls)
        if (iscell(tls))
            set(tls{t},'FontSize',fontsize);
        else 
            set(tls,'FontSize',fontsize);
        end
    end
end

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',dims);
set(gcf,'PaperPosition',[0 0 dims]);