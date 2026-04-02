if ~strcmp(test,'strat')
    title(str,txtopt{:})
    l = legend(lgnd,'location',loc);
    l.Visible = 'on';
    if exist('auxlines','var')
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot([tstart tend],auxlines([1 1],:),'-');
    end
    ax = gca;
    ax.XScale = xscale;
    ax.YScale = yscale;
    hold off
    if exist('yrange','var')
        ylim(yrange);
    end
    set(gca,txtopt{:})
end
if strcmp(task,'plotsave')
    if relaxation_flag
        print('-depsc2',['./save/' test '_' relax_meth '_' str_plt '_' num2str(idx) '.eps'])
    else
        print('-depsc2',['./save/' test '_0' str_plt '.eps'])
    end

end