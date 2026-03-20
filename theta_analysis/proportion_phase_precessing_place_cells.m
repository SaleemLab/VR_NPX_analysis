
function proportion_phase_precessing_place_cells

sessions = data_folders;
session_names = fieldnames(sessions);
PP =plotting_parameters;

c = 1;
for p = 1 : length(session_names)
    folders = sessions.(sprintf('%s',cell2mat(session_names(p))));
    for s = 1: length(folders)
        cd(cell2mat(folders(s)))
        
        load('extracted_phase_precession_absolute_location.mat')
        
        for t = 1 : length(TPP) % for each track
            
            sig_idx = [find(TPP(t).circ_lin_PVAL_dir1 < 0.05) find(TPP(t).circ_lin_PVAL_dir2 < 0.05)];
            proportion = length(unique(sig_idx))/length(TPP(t).cell_id);
            protocol.(sprintf('%s','T',num2str(t)))(c) = proportion;
            protocol.(sprintf('%s','scores_T',num2str(t))){c} = [TPP(t).circ_lin_PVAL_dir1(TPP(t).circ_lin_PVAL_dir1 < 0.05) TPP(t).circ_lin_PVAL_dir2(TPP(t).circ_lin_PVAL_dir2 < 0.05)];
        end
        c =c +1;
    end
end

%%% Proportion of phase precessing cells
mat= nan(20,8);
c=20;
for i = 1 : 5
    mat(1:4,i) = protocol.T2(c-3:c);
    c = c-4;
end
mat(:,6) =protocol.T1';
mat(:,7) =protocol.T3';
mat(:,8) =protocol.T4';

col = [PP.T2(5,:);PP.T2(4,:);PP.T2(3,:);PP.T2(2,:);PP.T2(1,:);PP.T1;[0.3 0.3 0.3];[0.6 0.6 0.6]];
figure;
boxplot(mat,'PlotStyle','traditional','Color',col,'LabelOrientation','horizontal','Widths',0.5);%,'BoxStyle','filled'
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
tt = get(a,'tag');   % List the names of all the objects
idx = find(strcmpi(tt,'box')==1);  % Find Box objects
boxes =a(idx);
set(boxes,'LineWidth',2); % Set width
box off

for j = 1 : 8
    hold on
    plot(j,mat(:,j),'o','MarkerFaceColor','w','MarkerEdgeColor',[0.3 0.3 0.3])
end
c=1;

ylim([0 1])
ylabel('Proportion of phase precessing place cells')
xticks([]);
ax=gca;
ax.FontSize = 16;
f=gcf;
f.Name = 'Proportion of precessing place cells per protocol';


%%% Scores of phase precessing cells
mat3 = [];
mat3 = [[protocol.scores_T2{:}]'; [protocol.scores_T1{:}]'; [protocol.scores_T3{:}]'; [protocol.scores_T4{:}]'];
grp = [];
c=20;
for i = 1 : 5
    grp = [grp; ones(length([protocol.scores_T2{c-3:c}]),1)*i];
    c = c-4;
end
grp = [grp; ones(length([protocol.scores_T1{:}]),1)*6; ones(length([protocol.scores_T3{:}]),1)*7];
c=20;
for i = 1 : 5
    grp = [grp; ones(length([protocol.scores_T4{c-3:c}]),1)*(i+7)];
    c = c-4;
end

figure
beeswarm(grp,mat3,'sort_style','nosort','colormap',[flipud(PP.T2);PP.T1;PP.T1;flipud(PP.T2)],'dot_size',2,'overlay_style','ci','corral_style','rand');


mat2= nan(1484,12);
c=20;
for i = 1 : 5
    mat2(1:length([protocol.scores_T2{c-3:c}]),i) = [protocol.scores_T2{c-3:c}];
    c = c-4;
end
mat2(1:length([protocol.scores_T1{:}]),6) =[protocol.scores_T1{:}]';
mat2(1:length([protocol.scores_T3{:}]),7) =[protocol.scores_T3{:}]';
%mat2(:,8) =protocol.scores_T4';
c=20;
for i = 1 : 5
    mat2(1:length([protocol.scores_T4{c-3:c}]),i+7) = [protocol.scores_T4{c-3:c}];
    c = c-4;
end

PP =plotting_parameters;
col = [flipud(PP.T2);PP.T1;PP.T1;flipud(PP.T2)];
figure;
boxplot(mat2,'PlotStyle','traditional','Color',col,'LabelOrientation','horizontal','Widths',0.5);%,'BoxStyle','filled'
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
tt = get(a,'tag');   % List the names of all the objects
idx = find(strcmpi(tt,'box')==1);  % Find Box objects
boxes =a(idx);
set(boxes,'LineWidth',2); % Set width
box off
for j = 1 : 12
    hold on
    plot(j,mat2(:,j),'o','MarkerFaceColor','w','MarkerEdgeColor',[0.3 0.3 0.3])
end
c=1;

ylim([0 1])
ylabel('Proportion of phase precessing place cells')
xticks([]);
ax=gca;
ax.FontSize = 16;
f=gcf;
f.Name = 'Proportion of precessing place cells per protocol';




save_all_figures('X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Chapter 3\Phase_precession',[])






end
