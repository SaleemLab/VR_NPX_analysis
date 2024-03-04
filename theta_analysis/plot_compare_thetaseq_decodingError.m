% DECODING ERROR COMPARISON BETWEEN PROTOCOLS

function plot_compare_thetaseq_decodingError(bayesian_control)

PP = plotting_parameters;

sessions = data_folders;
session_names = fieldnames(sessions);
if ~isempty(bayesian_control)
    path1 = '\Bayesian controls\Only first exposure';
    path2 = '\Bayesian controls\Only re-exposure';
    sessions_1 = arrayfun(@(y) arrayfun(@(x) append(sessions.(sprintf('%s',session_names{y})){x},path1),...
        1:length(sessions.(sprintf('%s',session_names{y}))),'UniformOutput',0),1:length(session_names),'UniformOutput',0);
    sessions_2 = arrayfun(@(y) arrayfun(@(x) append(sessions.(sprintf('%s',session_names{y})){x},path2),...
        1:length(sessions.(sprintf('%s',session_names{y}))),'UniformOutput',0),1:length(session_names),'UniformOutput',0);
end


c =1;
t1_DecError=[];
t3_DecError=[];
t1_lap_DecError = [];
t2_lap_DecError = [];
for p = 1 : length(session_names)
    folders = sessions.(sprintf('%s',session_names{p}));
    t2_DecError{p} = [];
    t4_DecError{p} = [];
    
    for s = 1: length(folders)
        cd(folders{s})
        disp(folders{s})
        if ~isempty(bayesian_control)
            folders_F = cell2mat(sessions_1{1,p}(s));
            load([folders_F '\Theta\thetaSequences_decodingError.mat']);
        else
            load 'Theta\thetaSequences_decodingError.mat'
        end
        t1_DecError =  [t1_DecError mean([thetaSequences_decodingError(1).bayesian_decodingError(1).thetaSequence(:).median_DecodingError ...
            thetaSequences_decodingError(2).bayesian_decodingError(1).thetaSequence(:).median_DecodingError])];
        t2_DecError{p} =  [t2_DecError{p} mean([thetaSequences_decodingError(1).bayesian_decodingError(2).thetaSequence(:).median_DecodingError ...
            thetaSequences_decodingError(2).bayesian_decodingError(2).thetaSequence(:).median_DecodingError])];

        if ~isempty(bayesian_control)
            folders_R = cell2mat(sessions_2{1,p}(s));
            load([folders_R '\Theta\thetaSequences_decodingError.mat']);
        end
        t3_DecError =  [t3_DecError mean([thetaSequences_decodingError(1).bayesian_decodingError(3).thetaSequence(:).median_DecodingError ...
            thetaSequences_decodingError(2).bayesian_decodingError(3).thetaSequence(:).median_DecodingError])];

        t4_DecError{p} = [t4_DecError{p} mean([thetaSequences_decodingError(1).bayesian_decodingError(4).thetaSequence(:).median_DecodingError ...
            thetaSequences_decodingError(2).bayesian_decodingError(4).thetaSequence(:).median_DecodingError])];
    end
    
end

if ~isempty(bayesian_control)
    save('X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Chapter 3\Controls\Bayesian controls\All_ThetaSequences_decodingError.mat','t1_DecError','t2_DecError',...
        't2_DecError','t4_DecError','-v7.3')
else
    save('X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Chapter 3\Controls\All_ThetaSequences_decodingError.mat','t1_DecError','t2_DecError',...
        't2_DecError','t4_DecError','-v7.3')
end
% Concatenate
temp_mat = nan(19,13);
temp_mat(:,1) = t1_DecError;
temp_mat(:,8) = t3_DecError;
ct = 9:13;
for ii = 1 : length(t2_DecError)
    temp_mat(1:length(t2_DecError{ii}),ii+1) = t2_DecError{ii};
    temp_mat(1:length(t4_DecError{ii}),ct(ii)) = t4_DecError{ii};
end


 col = [PP.T1;PP.T2(1,:);PP.T2(2,:);PP.T2(3,:);PP.T2(4,:);PP.T2(5,:);[1 1 1];PP.T1;PP.T2(1,:);PP.T2(2,:);PP.T2(3,:);PP.T2(4,:);PP.T2(5,:)];
 xlabels = {'T1','T2-8','T2-4','T2-3','T2-2','T2-1',[],'R-T1','RT2-8','RT2-4','RT2-3','RT2-2','RT2-1'};
 
 figure 
 boxplot(temp_mat ,'PlotStyle','traditional','Color',col,'LabelOrientation','horizontal','Widths',0.5);%,'BoxStyle','filled'
 a = get(get(gca,'children'),'children');   % Get the handles of all the objects
 tt = get(a,'tag');   % List the names of all the objects
 idx1 = find(strcmp(tt,'Outliers'));
 delete(a(idx1))
 idx = find(strcmpi(tt,'box')==1);  % Find Box objects
 boxes = a(idx([7:13]));  % Get the children you need (boxes for first exposure)
 boxes2 = a(idx([1:6])); % Get the children you need (boxes for second exposure)
 set(boxes,'LineWidth',2); % Set width
 set(boxes2,'LineStyle',':'); % Set line style for re-exposure plots
 set(boxes2,'LineWidth',2); % Set width
 idx2 = [find(strcmp(tt,'Upper Whisker')) find(strcmp(tt,'Lower Whisker'))];
 set(a(idx2),'LineStyle','-'); % Set width
 set(a(idx2),'LineWidth',0.5); % Set width
 idx3 = [find(strcmp(tt,'Upper Adjacent Value')) find(strcmp(tt,'Lower Adjacent Value'))];
 set(a(idx3),'LineWidth',0.5); % Set width
 
 box off
 ylabel('Median Dec Error','FontSize',14)
 a = get(gca,'XTickLabel');
 set(gca,'XTickLabel',a,'fontsize',14)
 hold on
 for ii = 1 : size(temp_mat,2)
        plot(ii,temp_mat(:,ii),'o','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',3)
 end 
 
[p_scores,tbl_scores,stats_scores] = kruskalwallis(temp_mat(:,1:6));
c = multcompare(stats_scores,'ctype','dunn-sidak','Display','off');  % if anova pval is < 0.05, run multiple comparisons
 
[p_scores,stats_scores] = ranksum(temp_mat(:,1),temp_mat(:,6));
end