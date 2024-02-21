function save_all_figures(save_path,filename)
% Save figures in .fig and .png format

figlist = findobj('type','figure');

for i = 1 : numel(figlist)
    if exist(filename,'var')
        name = filename;
    else 
        name = get(figlist(i),'Name');
        if contains(name,'/') 
            name(strfind(name,'/')) = ';';
        elseif contains(name,'\')
            name(strfind(name,'\')) = ';';
        elseif contains(name,':')
            name(strfind(name,':')) = '-';
        end
        if isempty(name)
           disp('Figure has no name!')
           name = ['Figure_' num2str(i)];
           if exist(name,'file')
                  name = ['Figure_' num2str(i) '_A'];
           end
        end
    end 
    if exist([save_path,'\figures'], 'dir')
        saveas(figlist(i),[save_path,'\figures\',name,'.fig']);
%         saveas(figlist(i),[save_path,'\figures\png_figs\',name,'.png']);
        saveas(figlist(i),[save_path,'\figures\png_figs\',name,'.pdf']);
        exportgraphics(figlist(i),[save_path,'\figures\png_figs\',name,'.pdf'],'ContentType','vector')
    else 
        disp('not saved in figures folder!')
        saveas(figlist(i),[save_path,'\',name,'.fig']);
%         saveas(figlist(i),[save_path,'\',name,'.png']);
        saveas(figlist(i),[save_path,'\',name,'.pdf']);
        exportgraphics(figlist(i),[save_path,'\',name,'.pdf'],'ContentType','vector')
    end
    close  
end

end