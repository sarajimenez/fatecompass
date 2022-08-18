%{
Transcription factor activity estimation using a linear model of gene regulation
and a novel framework for regularizing using data diffusion

Data: Veres et. al 2019 
System: Differentiation of hESCs to beta-like cells

Author: Sara Jimenez
Email: sarajcorrea@gmail.com
Molina's lab and Gradwohl's lab
IGBMC | Strasbourg University
%}

clc; clear; close all; 

% CAREFUL: this must be an absolute path --> change accordingly
fcnpath = '/Users/jimenezs/Documents/Doctorate/fatecompass/fatecompass/functions';
userpath(fcnpath)

%% Import pre-processed data 

datapath = '../../../data/human/';

% Reading data
umap2  = importdata(strcat(datapath,'umap_2.txt'));
umap10 = importdata(strcat(datapath,'umap_10.txt'));
G      = importdata(strcat(datapath,'normalized_gene_expression.tsv'));
G      = G.data'; % [cellsxgenes]
G      = sparse(G);
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
N  = N - Nm;

%% Selecting training (tr) and test (ts) subsets 

genestr = rand(numgenes,1)<0.8;
genests = ~genestr;

Ntr = N(genestr,:);
Nts = N(genests,:);
Etr = G(:,genestr);
Ets = G(:,genests);
numgenestr = sum(genestr);

%% Nearest neighbor graph for regularization using data diffusion

numberconnections = 61;
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

% plot(digraph(M),'XData',umap2(:,1),'YData',umap2(:,2))
% axis off
% 
% filename = strcat({'../../../output/human/figures/graph_data_diffusion_regularization.png'});
% saveas(gcf, char(filename))

M = sparse(M);

%% Estimating TF activities using minimum norm least-squares 

Atr = lsqminnorm(Ntr,Etr'); % [factorsxcells] 
Atr = sparse(Atr);

%% Cross-validation scheme to fit the value of t that minimizes the MSE

% This is run only once to get the value of t

Ats = Atr'; % [cellsxfactors]
maxt = 15;
MSE_tr = zeros(maxt+1,1);
MSE_ts = zeros(maxt+1,1);
MSE_tr(1) = mean(mean((Etr' - Ntr*Ats').^2));
MSE_ts(1) = mean(mean((Ets' - Nts*Ats').^2));

for t = 1:maxt
    Ats = M * Ats;
    MSE_tr(t+1) = mean(mean((Etr' - Ntr*Ats').^2));
    MSE_ts(t+1) = mean(mean((Ets' - Nts*Ats').^2)); 
    fprintf('%i %.5f %.5f\n',t,MSE_tr(t),MSE_ts(t))
    
end

figure()
plot(zscore(MSE_tr),'LineWidth',2); hold on
plot(zscore(MSE_ts),'LineWidth',2)
legend('Training','test')
set(gca,'FontSize',15)

filename = strcat({'../../../output/human/figures/zscore_MSEvst.png'});
saveas(gcf, char(filename))

figure()
plot(MSE_tr,'LineWidth',2); hold on
plot(MSE_ts,'LineWidth',2)
xline(5,'--r','LineWidth',1,'Label','Diffusion time 3','FontSize',15,'Color','k')
legend('Training','Test')
set(gca,'FontSize',15)

filename = strcat({'../../../output/human/figures/MSEvst.png'});
saveas(gcf, char(filename))

%% Regularizing activities using optimum t = 3 

Areg = M^3 * Atr'; % [cellsxfactors]

figure()
histogram(Areg,339);
set(gca,'FontSize',15)

filename = strcat({'../../../output/human/figures/tf_activities_distribution.png'});
saveas(gcf, char(filename))

%% Export activities to plot them over differentiation trajectories
% This file is read by the script differentiation_trajectories.m

Areg = full(Areg);
Areg2export = table(strrep(motifs, '"',''), Areg'); % [factorsxcells] 

writetable(Areg2export,'../../../output/human/regularized_activities.csv',...
    'WriteVariableNames',0);

%%
