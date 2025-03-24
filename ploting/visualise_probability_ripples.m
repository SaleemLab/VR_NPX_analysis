
   nprobe = 1
   [~,x,~] = calculate_event_probability(min(ripples(2).sharp_wave_peaktimes)',slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
figure;hold on;plot(mean(x))

%    [~,x,~] = calculate_event_probability(min(ripples(2).sharp_wave_peaktimes)',temp.ints.UP(temp.ints.UP(:,2)-temp.ints.UP(:,1)<2),-0.5:0.02:0.5,0);
% plot(mean(x))
% % % 
   nprobe = 2
   [~,x,~] = calculate_event_probability(min(ripples(2).sharp_wave_peaktimes)',slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
hold on;plot(mean(x))


   nprobe = 1
   [~,x,~] = calculate_event_probability(min(ripples(1).sharp_wave_peaktimes)',slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
figure;hold on;plot(mean(x))

%    [~,x,~] = calculate_event_probability(min(ripples(1).sharp_wave_peaktimes)',temp.ints.UP(temp.ints.UP(:,2)-temp.ints.UP(:,1)<2),-0.5:0.02:0.5,0);
% plot(mean(x))
% % % 
   nprobe = 2
   [~,x,~] = calculate_event_probability(min(ripples(1).sharp_wave_peaktimes)',slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
hold on;plot(mean(x))
nprobe = 1




[~,x,~] = calculate_event_probability(ripples(1).SWS_peaktimes,slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
figure;hold on;plot(mean(x))

[~,x,~] = calculate_event_probability(ripples(1).SWS_peaktimes,temp.ints.UP(temp.ints.UP(:,2)-temp.ints.UP(:,1)<2),-0.5:0.02:0.5,0);
plot(mean(x))
%
nprobe = 2
[~,x,~] = calculate_event_probability(ripples(1).SWS_peaktimes,slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
hold on;plot(mean(x))




nprobe = 1
[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
figure;hold on;plot(mean(x))

[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,temp.ints.UP(temp.ints.UP(:,2)-temp.ints.UP(:,1)<2),-0.5:0.02:0.5,0);
plot(mean(x))
%
nprobe = 2
[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,slow_waves(nprobe).UP_ints(slow_waves(nprobe).UP_ints(:,2)-slow_waves(nprobe).UP_ints(:,1)<2),-0.5:0.02:0.5,0);
hold on;plot(mean(x))





nprobe = 1
[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,slow_waves(nprobe).DOWN_ints,-0.5:0.02:0.5,0);
figure;hold on;plot(mean(x))

[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,temp.ints.DOWN,-0.5:0.02:0.5,0);
plot(mean(x))
%
nprobe = 2
[~,x,~] = calculate_event_probability(ripples(2).SWS_peaktimes,slow_waves(nprobe).DOWN_ints,-0.5:0.02:0.5,0);
hold on;plot(mean(x))