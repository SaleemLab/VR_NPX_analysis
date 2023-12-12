% Batch function to generate data structures from collated files for Temporal Frequency Flickering Stimulus %
function BatchLFPData_TF

%Define Data Paths
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
    ROOTPATH = 'X:\ibn-vision';
end

Datapath = fullfile(ROOTPATH,'DATA','SUBJECTS');

%State Subject(Mouse No) List
SubjectList = {'M22009','M22021','M22023','M22020'};


for i = 1:size(SubjectList,2)
    
    DataStruct = struct();
    
    thisAnimal = SubjectList{1,i}
    
    DataStruct.Mouse = thisAnimal;
    
    %File Directories
    MatPath = fullfile(Datapath,thisAnimal,'processed',[thisAnimal '.mat']);   
    TempDir = dir(fullfile(Datapath,thisAnimal,'processed'));
    Session = TempDir(3).name;
    StimInfoPath = dir(fullfile(Datapath,thisAnimal,'processed',TempDir(3).name,'*.csv'));
    DataInfo = readtable(fullfile(StimInfoPath.folder,StimInfoPath.name),'HeaderLines',0,'ReadVariableNames', true);
    Ephys_Datapath = fullfile(Datapath,thisAnimal,'ephys',Session);

    GNo = DataInfo.gNumber(ismember(DataInfo.StimulusName,'StaticTF'));
        
    folder = dir(fullfile(Ephys_Datapath, ['*','g', sprintf('%d',GNo)]));
    GLXDir = fullfile( folder.folder, folder.name);
    clear folder
    folder = dir(fullfile(GLXDir, [thisAnimal, '*']));
    GLXDir = fullfile(folder.folder, folder.name);
    clear folder
    
    [ChannelMapData] = SGLXMetaToCoords_ChannelMap_FR(GLXDir);

    %Load the data 
    DataTemp = load(MatPath);          
        
    TrialInfo = DataTemp.tstimData.StimIndex;
    Trial_1 = TrialInfo==1;
    Trial_2 = TrialInfo==2;
    Trial_4 = TrialInfo==4;
    Trial_8 = TrialInfo==8;
    
    StimTime = (DataTemp.ttimeVector>=0) & (DataTemp.ttimeVector<=3);
    PostStimTime = DataTemp.ttimeVector>3;
    
    StimTimeVector = DataTemp.ttimeVector(StimTime);
    PostStimTimeVector = DataTemp.ttimeVector(PostStimTime);
    
    % Get data arrays 
    
    StimResps_1 = DataTemp.resps(:,StimTime,Trial_1);
    StimResps_2 = DataTemp.resps(:,StimTime,Trial_2);
    StimResps_4 = DataTemp.resps(:,StimTime,Trial_4);
    StimResps_8 = DataTemp.resps(:,StimTime,Trial_8);
    
    PostStimResps_1 = DataTemp.resps(:,PostStimTime,Trial_1);
    PostStimResps_2 = DataTemp.resps(:,PostStimTime,Trial_2);
    PostStimResps_4 = DataTemp.resps(:,PostStimTime,Trial_4);
    PostStimResps_8 = DataTemp.resps(:,PostStimTime,Trial_8);  
    
    DataStruct.StimResps_1 = StimResps_1;
    DataStruct.StimResps_2 = StimResps_2;
    DataStruct.StimResps_4 = StimResps_4;
    DataStruct.StimResps_8 = StimResps_8;
    
    DataStruct.PostStimResps_1 = PostStimResps_1;
    DataStruct.PostStimResps_2 = PostStimResps_2;
    DataStruct.PostStimResps_4 = PostStimResps_4;
    DataStruct.PostStimResps_8 = PostStimResps_8;
    
    % Organise data arrays on the basis of the channel mapping
    % (Thalamus - V1)
    
    Block1 = ChannelMapData.ycoords<1000;
    Block2 = ChannelMapData.ycoords>1000;
    
    [~,ShankOrderVT] = sort(ChannelMapData.shankind(Block1),'ascend');
    [~,ShankOrderV1] = sort(ChannelMapData.shankind(Block2),'ascend');
           
    StimRespsVT_1 = StimResps_1(Block1,:,:);
    StimRespsVT_2 = StimResps_2(Block1,:,:);
    StimRespsVT_4 = StimResps_4(Block1,:,:);
    StimRespsVT_8 = StimResps_8(Block1,:,:);
    
    StimRespsVT_1 = StimRespsVT_1(ShankOrderVT,:,:);
    StimRespsVT_2 = StimRespsVT_2(ShankOrderVT,:,:);
    StimRespsVT_4 = StimRespsVT_4(ShankOrderVT,:,:);
    StimRespsVT_8 = StimRespsVT_8(ShankOrderVT,:,:);
    
    StimRespsV1_1 = StimResps_1(Block2,:,:);
    StimRespsV1_2 = StimResps_2(Block2,:,:);
    StimRespsV1_4 = StimResps_4(Block2,:,:);
    StimRespsV1_8 = StimResps_8(Block2,:,:);
    
    StimRespsV1_1 = StimRespsV1_1(ShankOrderV1,:,:);
    StimRespsV1_2 = StimRespsV1_2(ShankOrderV1,:,:);
    StimRespsV1_4 = StimRespsV1_4(ShankOrderV1,:,:);
    StimRespsV1_8 = StimRespsV1_8(ShankOrderV1,:,:);
    
    PostStimRespsVT_1 = PostStimResps_1(Block1,:,:);
    PostStimRespsVT_2 = PostStimResps_2(Block1,:,:);
    PostStimRespsVT_4 = PostStimResps_4(Block1,:,:);
    PostStimRespsVT_8 = PostStimResps_8(Block1,:,:);
    
    PostStimRespsVT_1 = PostStimRespsVT_1(ShankOrderVT,:,:);
    PostStimRespsVT_2 = PostStimRespsVT_2(ShankOrderVT,:,:);
    PostStimRespsVT_4 = PostStimRespsVT_4(ShankOrderVT,:,:);
    PostStimRespsVT_8 = PostStimRespsVT_8(ShankOrderVT,:,:);
    
    PostStimRespsV1_1 = PostStimResps_1(Block2,:,:);
    PostStimRespsV1_2 = PostStimResps_2(Block2,:,:);
    PostStimRespsV1_4 = PostStimResps_4(Block2,:,:);
    PostStimRespsV1_8 = PostStimResps_8(Block2,:,:);
    
    PostStimRespsV1_1 = PostStimRespsV1_1(ShankOrderV1,:,:);
    PostStimRespsV1_2 = PostStimRespsV1_2(ShankOrderV1,:,:);
    PostStimRespsV1_4 = PostStimRespsV1_4(ShankOrderV1,:,:);
    PostStimRespsV1_8 = PostStimRespsV1_8(ShankOrderV1,:,:);
    
    figure
    title(['Thalamus responses during stimulus ', thisAnimal])
    subplot(221), imagesc(StimTimeVector, 1:384, mean(StimRespsVT_1,3))
    subplot(222), imagesc(StimTimeVector, 1:384, mean(StimRespsVT_2,3))
    subplot(223), imagesc(StimTimeVector, 1:384, mean(StimRespsVT_4,3))
    subplot(224), imagesc(StimTimeVector, 1:384, mean(StimRespsVT_8,3))
      
    figure
    title(['Thalamus responses post stimulus ', thisAnimal])
    subplot(221), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsVT_1,3))
    subplot(222), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsVT_2,3))
    subplot(223), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsVT_4,3))
    subplot(224), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsVT_8,3))
    
    figure
    title(['V1 responses during stimulus ', thisAnimal])
    subplot(221), imagesc(StimTimeVector, 1:384, mean(StimRespsV1_1,3))
    subplot(222), imagesc(StimTimeVector, 1:384, mean(StimRespsV1_2,3))
    subplot(223), imagesc(StimTimeVector, 1:384, mean(StimRespsV1_4,3))
    subplot(224), imagesc(StimTimeVector, 1:384, mean(StimRespsV1_8,3))
      
    figure
    title(['V1 responses post stimulus ', thisAnimal])
    subplot(221), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsV1_1,3))
    subplot(222), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsV1_2,3))
    subplot(223), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsV1_4,3))
    subplot(224), imagesc(PostStimTimeVector, 1:384, mean(PostStimRespsV1_8,3))
    
    figure
    title(['Average responses across all chanels ', thisAnimal])
    subplot(421)
    hold on
    plot(StimTimeVector,mean(mean(StimRespsVT_1,3),1))
    plot(StimTimeVector,mean(mean(StimRespsV1_1,3),1))
   
    subplot(422)
    hold on
    plot(PostStimTimeVector,mean(mean(PostStimRespsVT_1,3),1))
    plot(PostStimTimeVector,mean(mean(PostStimRespsV1_1,3),1))
    
    subplot(423)
    hold on
    plot(StimTimeVector,mean(mean(StimRespsVT_2,3),1))
    plot(StimTimeVector,mean(mean(StimRespsV1_2,3),1))
   
    subplot(424)
    hold on
    plot(PostStimTimeVector,mean(mean(PostStimRespsVT_2,3),1))
    plot(PostStimTimeVector,mean(mean(PostStimRespsV1_2,3),1))
    
    subplot(425)
    hold on
    plot(StimTimeVector,mean(mean(StimRespsVT_4,3),1))
    plot(StimTimeVector,mean(mean(StimRespsV1_4,3),1))
   
    subplot(426)
    hold on
    plot(PostStimTimeVector,mean(mean(PostStimRespsVT_4,3),1))
    plot(PostStimTimeVector,mean(mean(PostStimRespsV1_4,3),1))
    
    subplot(427)
    hold on
    plot(StimTimeVector,mean(mean(StimRespsVT_8,3),1))
    plot(StimTimeVector,mean(mean(StimRespsV1_8,3),1))
   
    subplot(428)
    hold on
    plot(PostStimTimeVector,mean(mean(PostStimRespsVT_8,3),1))
    plot(PostStimTimeVector,mean(mean(PostStimRespsV1_8,3),1))
    
    DataStruct.VT_Stim{1,1} = StimRespsVT_1;
    DataStruct.VT_Stim{2,1} = StimRespsVT_2;
    DataStruct.VT_Stim{3,1} = StimRespsVT_4;
    DataStruct.VT_Stim{4,1} = StimRespsVT_8;
    
    DataStruct.VT_PostStim{1,1} = PostStimRespsVT_1;
    DataStruct.VT_PostStim{2,1} = PostStimRespsVT_2;
    DataStruct.VT_PostStim{3,1} = PostStimRespsVT_4;
    DataStruct.VT_PostStim{4,1} = PostStimRespsVT_8;
    
    DataStruct.V1_Stim{1,1} = StimRespsV1_1;
    DataStruct.V1_Stim{2,1} = StimRespsV1_2;
    DataStruct.V1_Stim{3,1} = StimRespsV1_4;
    DataStruct.V1_Stim{4,1} = StimRespsV1_8;
    
    DataStruct.V1_PostStim{1,1} = PostStimRespsV1_1;
    DataStruct.V1_PostStim{2,1} = PostStimRespsV1_2;
    DataStruct.V1_PostStim{3,1} = PostStimRespsV1_4;
    DataStruct.V1_PostStim{4,1} = PostStimRespsV1_8;
    
    DataStruct.StimTimeVector = StimTimeVector;
    DataStruct.PostStimTimeVector = PostStimTimeVector;
    
    save(['X:\ibn-vision\DATA\PROJECTS\Npx_TemporalFrequency\', thisAnimal], 'DataStruct', '-v7.3');
          
    clearvars -except ROOTPATH Datapath SubjectList 
end

