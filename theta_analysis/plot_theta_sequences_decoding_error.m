% PLOT DECODING ERROR IN THETA WINDOWS
% MH, 2020
% Plots decoding error per track (rows) and directions (columns)

function plot_theta_sequences_decoding_error(save_option)


load Theta\thetaSequences_decodingError.mat

% For visualization purposes, create a structure concatenating 3 cycles together

for d = 1 : length(thetaSequences_decodingError) % for each direction
    for t = 1 : length(thetaSequences_decodingError(d).bayesian_decodingError) % for each track
        a = {thetaSequences_decodingError(d).bayesian_decodingError(t).thetaSequence(:).decodingErrors};
        A = a(bsxfun(@plus,(1:3),(0:1:length(a)-3)')); %concatenate 3 cycles together, overlapping every two (123,234,345..)
        for i = 1 : size(A,1)
            concat_DecErrors.(sprintf('%s','direction',num2str(d)))(t).track(i,:) = cat(2,A{i,:}); %merge together each 3 cycles in a row
        end
        clear a A 
    end
end

PP = plotting_parameters;

folder = pwd;
protocol = str2num(folder(end));
protocols = [8,4,3,2,1];
curr_prot = find(protocols == protocol);

name = strsplit(pwd,'\');
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
f1.Name = ['ThetaSeq_DecError_' name{end}];
c = [1,3,5,7,2,4,6,8];
i =1;

for d = 1 : length(thetaSequences_decodingError) % for each direction
    
    for t = 1 : length(thetaSequences_decodingError(d).bayesian_decodingError) % for each track
        average_decoding_errors = mean(concat_DecErrors.(sprintf('%s','direction',num2str(d)))(t).track,1);
        
        ax(c(i)) = subplot(4,2,c(i));
        hold on
        plot(average_decoding_errors,'Color',PP.P(curr_prot).colorT(t,:),'LineWidth',PP.Linewidth{t},'LineStyle',PP.Linestyle{t})
        xticks([2,7,12,17,22,27])
        xticklabels({'-2\pi','-\pi','0','\pi','2\pi','5\pi/2'})
        xlabel('Phase')
        ylabel('Decoding error (cm)') ;
        title(['Track ' num2str(t)])
        ax(c(i)).FontSize = 12;
        
        i = i+1;
    end
    
end

linkaxes([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6) ax(7) ax(8)],'y')

for jj = 1 : length(c)
    plot(ax(jj),[7 7],[min(ylim) max(ylim)],':','Color',[0.5 0.5 0.5],'LineWidth',2)
    plot(ax(jj),[17 17],[min(ylim) max(ylim)],':','Color',[0.5 0.5 0.5],'LineWidth',2)
    plot(ax(jj),[27 27],[min(ylim) max(ylim)],':','Color',[0.5 0.5 0.5],'LineWidth',2)    
end

if save_option == 1
    save_fol = 'X:\BendorLab\Drobo\Lab Members\Marta\Analysis\HIPP\Chapter 3\Controls';
    save_all_figures(save_fol,[])
end

end 