function    [tuning_curve] = plot_tuning_curve(grating_response,if_plot,if_plot_psth)


no_unit = size(grating_response{1},1);
orientation_angles = 0:15:165;
angle_color = hsv(12);
time_vec = 180/(12*150):180/(12*150):180;
time_vec = reshape(time_vec,[150 12]);

MUA_filter_length = 30;
SD_alpha = 2; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

mean_diff_resp = cell(12,no_unit);
avg_mean_diff_resp = zeros(12,no_unit);
se_mean_resp = zeros(12,no_unit);
for iFig = 1:ceil(no_unit/6)
    if if_plot
        figure;
    end
    for iUnit = 1:6
        if iUnit+iFig*6 > no_unit
            break
        end
        if if_plot
            subplot(2,3,iUnit)
        
        hold on
        end
        for iAngle = 1:length(orientation_angles)
            resp_tmp = squeeze(grating_response{iAngle}(iUnit+(iFig-1)*6,:,:));
            smooth_resp = filtfilt(w,1,resp_tmp);
            %                     plot(repmat(time_vec(:,iAngle),[1,size(resp_tmp,2)]),smooth_resp,'color', [.5 .5 .5]);
            %                     alpha(.5)
            %average response and standard deviation
            mean_resp = mean(smooth_resp,2);
            std_resp = std(smooth_resp,[],2);
            se_resp = std_resp./sqrt(size(resp_tmp,2));
            if if_plot_psth
                plot(time_vec(:,iAngle),mean(smooth_resp,2),'color',angle_color(iAngle,:));
                fill([time_vec(:,iAngle); flipud(time_vec(:,iAngle))], [(mean_resp-se_resp); flipud((mean_resp+se_resp))], angle_color(iAngle,:),'linestyle','none','FaceAlpha', 0.2,'HandleVisibility','off');
                line(time_vec(:,iAngle),mean_resp,'Color', angle_color(iAngle,:), 'LineWidth', 1.5);
            end
            %tuning curve taking from mean of the 1s stimulus onset
            %time
            mean_on_resp = mean(smooth_resp(31:90,:),1);
            mean_off_resp = mean(smooth_resp(91:end,:),1);
            mean_diff_resp{iAngle,iUnit+(iFig-1)*6} = abs(mean_on_resp - mean_off_resp);
            avg_mean_diff_resp(iAngle,iUnit+(iFig-1)*6)=mean(mean_diff_resp{iAngle,iUnit+(iFig-1)*6});
            se_mean_resp(iAngle,iUnit+(iFig-1)*6) = std(mean_diff_resp{iAngle,1})./sqrt(size(resp_tmp,2));
            if if_plot
                scatter(repmat((iAngle-1)*15,[length(mean_resp),1]),mean_diff_resp{iAngle,iUnit+(iFig-1)*6},[],angle_color(iAngle,:));
                title(num2str(iUnit+(iFig-1)*6))
            end
        end
        if if_plot
            errorbar(0:15:165,avg_mean_diff_resp,se_mean_resp);
            hold off
        end
    end
end

tuning_curve.mean_diff_resp = mean_diff_resp;
tuning_curve.avg_mean_diff_resp = avg_mean_diff_resp;
tuning_curve.se_resp = se_mean_resp;

end



