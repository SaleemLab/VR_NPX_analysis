% function AssemblyTemplates = assembly_patterns(SpikeCount,opts)
function [AssemblyTemplates,AssemblyIDs,time_projection,ASSEMBLYPROJECTOR,variance] = assembly_patterns(SpikeCount, varargin)

%assembly_patterns - extracts assembly patterns from the spike matrix.
%Patterns = assembly_patterns(Activitymatrix,opts)
% 
% INPUTS	
% SpikeCount      spike matrix. Rows represent neurons, columns represent 
%                 time bins. Thus, each element of the matrix carries the 
%                 spike count of a given neuron at a given bin.
%    <options>   optional list of property-value pairs (see table below)
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     opts          set parameters. All fields described below.
%       opts.threshold.method: defines the method to compute the threshold 
%                              for assembly detection. Options are:   
%         'MarcenkoPastur': uses the analytical bound.
%         'binshuffling': estimate eigenvalue distribution for independent 
%                         activity from surrogate matrices generated by 
%                         shuffling time bins.
%         'circularshift': estimate eigenvalue distribution for independent 
%                          activity from surrogate matrices generated by 
%                          random circular shifts of original spike matrix.
%       opts.threshold.permutations_percentile: defines which percentile of 
%                              the surrogate distribution of maximal 
%                              eigenvalues is used as statistical 
%                              threshold. It must be a number between 0 and 
%                              100 (95 or larger recommended). Not used 
%                              when 'MarcenkoPastur' is chosen.
%       opts.threshold.number_of_permutations: defines how many surrogate 
%                              matrices are generated (100 or more 
%                              recommended). Not used when 'MarcenkoPastur' 
%                              is chosen.
%       opts.Patterns.method: defines which method is used to extract 
%                             assembly patterns. Options are: 'PCA' or 
%                             'ICA' (recommended). 
%       opts.Patterns.number_of_iterations: number of iterations for 
%                             fastICA algorithm (100 or more recommended). 
%                             Not used when 'PCA' is chosen.
%    =========================================================================
%
% OUTPUTS
%    Patterns       assembly patterns. Columns denote assembly # and rows 
%                   neuron #.
% 
% EXAMPLES
% 
% Patterns = assembly_patterns(Activitymatrix); The following default 
% options will be used: 
%   opts.Patterns.method: 'ICA'
%   opts.threshold.method: 'MarcenkoPastur'
%   opts.Patterns.number_of_iterations: 500.
% 
% This framework is described in: 
% Lopes-dos-Santos V, Ribeiro S, Tort ABL 
% (2013) Detecting cell assemblies in large neuronal populations, 
% Journal of Neuroscience Methods.
%
% Adapted to neurocode by aza

% Parse inputs 
p=inputParser;
% addParameter(p,'SpikeCount',[],@isnumeric);
addParameter(p,'n_perm',100,@isnumeric);
addParameter(p,'n_perm_percentile',95,@isnumeric);
addParameter(p,'n_iter',500,@isnumeric);
addParameter(p,'n_SD',2,@isnumeric);
addParameter(p,'patterns_method','ICA',@ischar);
addParameter(p,'threshold_method','MarcenkoPastur',@ischar);

parse(p,varargin{:})
% SpikeCount = p.SpikeCount;
n_perm = p.Results.n_perm;
n_perm_percentile = p.Results.n_perm_percentile;
n_iter = p.Results.n_iter;
n_SD = p.Results.n_SD;
patterns_method = p.Results.patterns_method;
threshold_method = p.Results.threshold_method;

%z-score spike count
zSpikeCount = zscore(SpikeCount');

%correlation matrix
CorrMatrix = corr(zSpikeCount);
CorrMatrix(isnan(CorrMatrix))=0;
[eigenvectors,d] = eig(CorrMatrix);
eigenvalues=diag(d);

q = size(zSpikeCount,1)/size(zSpikeCount,2);

if q<1
    error('Error: Number of time bins must be larger than number of neurons \n');
end

%% Finding number of assemblies

switch threshold_method
    case 'MarcenkoPastur'
        fprintf('Using Marcenko-Pastur distribution for estimating number of assemblies \n')
        lambda_max = ((1+sqrt(1/q))^2);
    case 'binshuffling'
        fprintf('Generating control spike count matrix for estimating number of assemblies \n')
%         if ~isfield(opts.threshold,'number_of_permutations')
%             auxmsdg = 'Please enter a number of surrogates larger than zero: ';
%             opts.threshold.number_of_permutations = input(auxmsdg);
%         end
%         if ~isfield(opts.threshold,'permutations_percentile')
%             auxmsdg = 'Please enter percentile for statistical threshold: ';
%             opts.threshold.permutations_percentile = input(auxmsdg);
%         end
%         while opts.threshold.number_of_permutations<=0
%             auxmsdg = 'Please enter a number of surrogates larger than zero: ';
%             opts.threshold.number_of_permutations = input(auxmsdg);
%         end
        fprintf(['Number of permutations:  ' num2str(n_perm) '\n'])
        control_max_eig = bin_shuffling(SpikeCount,n_perm);
        lambda_max = prctile(control_max_eig,n_perm_percentile);
    case 'circularshift'
        fprintf('Generating control spike count matrix for estimating number of assemblies \n')
%         if ~isfield(opts.threshold,'number_of_permutations')
%             auxmsdg = 'Please enter a number of surrogates larger than zero: ';
%             opts.threshold.number_of_permutations = input(auxmsdg);
%         end
%         while opts.threshold.number_of_permutations<=0
%             auxmsdg = 'Please enter a number of surrogates larger than zero: ';
%             opts.threshold.number_of_permutations = input(auxmsdg);
%         end
%         if ~isfield(opts.threshold,'permutations_percentile')
%             auxmsdg = 'Please enter percentile for statistical threshold: ';
%             opts.threshold.permutations_percentile = input(auxmsdg);
%         end
        fprintf(['Number of permutations:  ' num2str(n_perm) '\n'])
        control_max_eig = circular_shift(SpikeCount,n_perm);
        lambda_max = prctile(control_max_eig,n_perm_percentile);
end
NumberOfAssemblies = sum(eigenvalues>lambda_max);
fprintf(['Number of assemblies detected: ' num2str(NumberOfAssemblies) '\n'])
if NumberOfAssemblies<1
    AssemblyTemplates=[];
    returnnonsignifassemblies = eigenvectors(:,eigenvalues<=lambda_max);
end
%% Finding co-activation patterns

switch patterns_method
    case 'PCA'
        [garbage,PC_position] = sort(-eigenvalues);
        AssemblyTemplates = eigenvectors(:,PC_position(1:NumberOfAssemblies));

    case 'ICA'
        AssemblyTemplates=...
            fast_ica(zSpikeCount,NumberOfAssemblies,n_iter);
end


% The sign of the weights (AssemblyTemplates) in a component is arbitrary (+component and -component are equivalent)
% Flip weights so that the most deviating weight of a component is positive (it is more convenient for visualisation that the assembly has positive weights)
flip = max(AssemblyTemplates)<-min(AssemblyTemplates);
AssemblyTemplates(:,flip) = -AssemblyTemplates(:,flip);

% Normalise weights as Van de Ven et al (2016):
nor = ones(1,size(AssemblyTemplates,2)); for i=1:size(AssemblyTemplates,2),nor(1,i)=norm(AssemblyTemplates(:,i)); end
AssemblyTemplates = bsxfun(@rdivide,AssemblyTemplates,nor);


signifassemblies = eigenvectors(:,eigenvalues>lambda_max);
nonsignifassemblies = eigenvectors(:,eigenvalues<=lambda_max);
thresh = mean(eigenvectors(:))+n_SD*std(eigenvectors(:));
AssemblyIDs.signifassemblies=signifassemblies;
AssemblyIDs.nonsignifassemblies=nonsignifassemblies;
AssemblyIDs.thresh=thresh;

%% Finding Activity Patterns
zSpikeCount=zSpikeCount';
time_projection=zeros(size(AssemblyTemplates,2),size(zSpikeCount,2));

variance = var(zSpikeCount'*AssemblyTemplates)/size(zSpikeCount',2);
[~,order] = sort(-variance); % from highest to lowest
AssemblyTemplates = AssemblyTemplates(:,order);
variance = variance(:,order);

for assembly_idx = 1:size(AssemblyTemplates,2)
    
    % computing projector
    ASSEMBLYPROJECTOR(:,:,assembly_idx)=AssemblyTemplates(:,assembly_idx)*AssemblyTemplates(:,assembly_idx)';
    ASSEMBLYPROJECTOR(:,:,assembly_idx)=squeeze(ASSEMBLYPROJECTOR(:,:,assembly_idx))-diag(diag(squeeze(ASSEMBLYPROJECTOR(:,:,assembly_idx))));
    
    % computing activity time course
    for ntime=1:size(zSpikeCount,2)
        
        time_projection(assembly_idx,ntime)=(zSpikeCount(:,ntime)'*ASSEMBLYPROJECTOR(:,:,assembly_idx)*zSpikeCount(:,ntime));
        
    end
    
end
