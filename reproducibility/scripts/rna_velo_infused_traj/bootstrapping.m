%{
Bootstrapping to get the empirical distribution of the 
transcription factor activities

Data: Bastidas-Ponce et. al 2019 
System: Mouse in-vivo at E15.5

Author: Sara Jimenez
Email: sarajcorrea@gmail.com
Molina's lab and Gradwohl's lab
IGBMC | Strasbourg University

%}

clc; clear; close all; 

% CAREFUL: this must be an absolute path --> change accordingly
fcnpath = '/Users/jimenezs/Documents/Doctorate/fatecompass/functions';
userpath(fcnpath)

%% Import pre-processed data 

datapath = '../../../data/mouse/';

% Reading data
umap2  = importdata(strcat(datapath,'umap_2.txt'));
umap10 = importdata(strcat(datapath,'umap_10.txt'));
G      = importdata(strcat(datapath,'normalized_gene_expression.tsv'));
G      = G.data'; % [cellsxgenes]
N      = importdata(strcat(datapath,'binding_sites.csv'));
motifs = split(N.textdata(1,1),";");
N      = N.data; % [genesxfactors]
N      = sparse(N);

% Selecting dimensions
[numcells,numgenes] = size(G);
[~,numfac]          = size(N);

%% Normalization for activity estimation

% Cell- and gene- normalized log-expression values
G = G - (1/numgenes)*sum(G,2);
G = G - (1/numcells)*sum(G,1);

% Motif-normalized site-counts
Nm = (1/numgenes)*sum(N);
N = N - Nm;

%% Nearest neighbor graph for regularization using data diffusion

numberconnections = 11;
I = knnsearch(umap10,umap10,'K',numberconnections);
M = zeros(numcells);
for i = 1:numcells
    M(i,I(i,2:end)) = 1;
end

% We want a cell's own observed values to have the highest impact on the 
% imputation of its own values; therefore, our transition matrix M^* allows for 
% self-loops, and these are the most probable steps in the random walk.
M = M + 10*eye(size(M)); 
M = M./sum(M,2);
M = M';

M = sparse(M);

%% Resampling 

numsamples  = 100; % Number of times to resample
Est_dist    = cell(numsamples,1);  

for i = 1:numsamples
    
    % resampling
    genestr = rand(numgenes,1)<0.8;
    Ntr = N(genestr,:);
    Etr = G(:,genestr);
    
    % Compute the estimate of the activities
    Atr = lsqminnorm(Ntr,Etr'); % [factorsxcells]
    
    % Save the regularized value of the activities
    Est_dist{i} = M^5 * Atr'; % [cellsxfactors]
    
    fprintf('%i\n',i)
        
end

%% Formating the distributions
% One entry per factor and one matrix per factor containing the different
% realizations of activities for that given factor across cells. 

mydist = cell(numfac,1);

for f = 1:numfac 
    
    for i = 1:numsamples
        
        mydist{f}(:,i) = Est_dist{i}(:,f);
        
    end
      
end

%% Plot empirical distribution for a given factor and a given cell

myfac = 10;
mycell = 200;

figure()
subplot(1,2,1)
histogram(mydist{myfac}(mycell,:),1000);
title('Distribution of the sample mean')

subplot(1,2,2)
ksdensity(mydist{myfac}(mycell,:));
title('Distribution estimate of the sample mean')

%% Export distributions to estimate z-score
% This file is read by the script differential_motif_activity_analysis.m

filename = '../../../output/mouse/bootstrap_distribution_activities.mat';
save(filename,'mydist')

%%

