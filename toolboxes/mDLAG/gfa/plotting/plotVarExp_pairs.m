function plotVarExp_pairs(pairVarExp, pairs, numGroups, varargin)
%
% plotVarExp_pairs(pairVarExp, pairs, numGroups, ...)
%
% Description: Visualize the shared variance explained in each group by
%              their pairwise interaction with another group.
% 
% Arguments:
%
%     Required:
%
%     pairVarExp -- (numPairs x 2) array;
%                       pairVarExp(i,1): shared variance explained by
%                                        pairwise interaction, group 1 of
%                                        pair i
%                       pairVarExp(i,2): shared variance explained by
%                                        pairwise interaction, group 2 of
%                                        pair i
%     pairs      -- (numPairs x 2) array; pairs(i,:) gives the indexes of
%                   groups in pair i, and corresponds to the relevant rows
%                   of pairDims and pairVarExp.
%     numGroups  -- int; number of observed groups
%
%     Optional:
%
%     pairVarExp_sem -- (numPairs x 2) array; The standard error of each 
%                       element in pairVarExp, if pairVarExp is the mean  
%                       over different runs. (default: [])
%     groupNames     -- (1 x numGroups) cell array; List of strings 
%                       containing the name of each group. By default, 
%                       names will be '1', '2', '3', etc.
%
% Outputs:
%
%     none.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Sep 2022 -- Initial full revision.
%     22 Sep 2022 -- Added option to plot SEM.

pairVarExp_sem = [];
groupNames = cell(1,numGroups);
for groupIdx = 1:numGroups
    groupNames{groupIdx} = num2str(groupIdx); 
end
assignopts(who, varargin);

numPairs = size(pairs,1);

% Set up x-axis labels
xtcklbls = cell(numPairs,3);
for pairIdx = 1:numPairs
    % Total in group 1
    xtcklbls{pairIdx,1} = sprintf('Total, %s', groupNames{pairs(pairIdx,1)});
    % Shared between both groups
    xtcklbls{pairIdx,2} = sprintf('%s-%s', groupNames{pairs(pairIdx,1)}, groupNames{pairs(pairIdx,2)});
    % Total in group 2
    xtcklbls{pairIdx,3} = sprintf('Total, %s', groupNames{pairs(pairIdx,2)});
end

% Plot shared variances explained
figure;
for pairIdx = 1:numPairs
    subplot(1,numPairs,pairIdx);
    hold on;
    bar(pairVarExp(pairIdx,:), ...
        'facealpha', 0.3, ...
        'facecolor', 'k', ...
        'edgecolor', 'k', ...
        'linewidth', 1.5);
    if ~isempty(pairVarExp_sem)
        errorbar(1:length(pairVarExp(pairIdx,:)), pairVarExp(pairIdx,:), pairVarExp_sem(pairIdx,:), ...
             'k.', ...
             'linewidth', 1.5, ...
             'marker', 'none'); 
    end
    xlabel('Group');
    ylabel('Shared var. exp. by interaction');
    xticks(1:length(pairVarExp(pairIdx,:)));
    xticklabels({groupNames{pairs(pairIdx,1)}, groupNames{pairs(pairIdx,2)}});
    title(xtcklbls{pairIdx,2});
    ylim([0 1]);
end
