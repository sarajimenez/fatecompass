%{
Simulation of stochastic trajectories during beta-cell differentiation 

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

umap2  = importdata(strcat(datapath,'umap_2.txt'));
umap10 = importdata(strcat(datapath,'umap_10.txt'));
cellID = importdata(strcat(datapath,'cell_id.csv'));
GE     = importdata(strcat(datapath,'normalized_gene_expression.tsv'));
G      = GE.data'; % [cellsxgenes]

% Selecting dimmensions
[numcells,numgenes] = size(G); 

% Containers to map cell-types and gene names
clusters = {'prog_nkx61','neurog3_early','neurog3_mid', 'neurog3_late',...
    'fev_high_isl_low','sst_hhex','ph','scbeta', 'ec'};
mapcellID = containers.Map(clusters,1:length(clusters));
mapgeneID = containers.Map(GE.textdata,1:length(GE.textdata));

%% Figure parameters

Prog_c  = [102 195 166]/255;
Ngn3E_c = [243 194 123]/255;
Ngn3M_c = [238 133 50]/255;
Ngn3L_c = [155 89 51]/255;
B_c     = [165 225 127]/255;
EC_c    = [238 132 124]/255;
A_c     = [53 122 197]/255;
D_c     = [107 56 154]/255;
Fev_c   = [201 177 215]/255;

colcluster = strcmp(cellID,'prog_nkx61')*Prog_c + ...
    strcmp(cellID,'neurog3_early')*Ngn3E_c + ...
    strcmp(cellID,'neurog3_mid')*Ngn3M_c + ...
    strcmp(cellID,'neurog3_late')*Ngn3L_c + ...
    strcmp(cellID,'fev_high_isl_low')*Fev_c + ...
    strcmp(cellID,'sst_hhex')*D_c + ...
    strcmp(cellID,'ph')*A_c + ...
    strcmp(cellID,'scbeta')*B_c + ...
    strcmp(cellID,'ec')*EC_c;

%% Nearest neighbor graph representing phenotypic manifold

numnb = 51; 
K = knnsearch(umap10,umap10,'K',numnb); % This gives indices 
K = K(:,2:end);
numnb = numnb - 1;

% Building network
P = zeros(numcells);
for i = 1:numcells
    P(i,K(i,:)) = 1; % adjancency matrix
end

Gnn = digraph(P);
plot(Gnn,'XData',umap2(:,1),'YData',umap2(:,2))
axis off

filename = strcat({'../../../output/human/figures/nearest_neighbor_graph.png'});
saveas(gcf, char(filename))

%% Pre-defining sources and sinks 

root = 2383;%1104;
betaSink = [15326;15959;16260;16669;16906;17078;17248;17444;17566;17623;17639;...
    17913;17918;18026;18113;18153;18370;18433;18614;18875;18905;18946;...
    18978;19028;19072;19113;19128;19149;19178;19203;19260;19275;19285;...
    19296;19372;19373;19409;19469;19491;19557;19571;19590;19635;19647;...
    19673;19676;19713;19798;19827;19878];
alphaSink = [14791;15096;15098;15538;15708;15714;15735;15881;16008;16098;16111;...
    16132;16139;16159;16288;16552;16805;16844;16847;16967;16974;17037;...
    17102;17146;17149;17173;17180;17187;17193;17252;17306;17349;17368;...
    17407;17425;17456;17462;17489;17513;17521;17580;17629;17736;17767;...
    17773;17809;17811;17845;17854;17994];
ecSink  = [8233;10714;14068;14250;14382;14420;14501;14527;14529;14588;14678;...
    14680;14697;14747;14750;14794;14832;14840;14893;14904;15214;15266;...
    15338;15341;15458;15522;15597;15618;15642;15762;15844;15905;15989;...
    16189;16228;16235;16239;16243;16245;16372;16407;16502;16544;16709;...
    16753;16775;16834;16942;16977;16982];

figure()
hold on
scatter(umap2(:,1),umap2(:,2),5,[190 190 190]/255)
scatter(umap2(alphaSink,1),umap2(alphaSink,2),10,A_c,'filled')
scatter(umap2(betaSink,1),umap2(betaSink,2),10,B_c,'filled')
scatter(umap2(ecSink,1),umap2(ecSink,2),10,EC_c,'filled')
scatter(umap2(root,1),umap2(root,2),50,Prog_c,'filled')
axis off

filename = strcat({'../../../output/human/figures/sources_sinks.png'});
saveas(gcf, char(filename))

%% Building energy landscape 

nodes = [root betaSink' alphaSink' ecSink'];
ntype = 1/length(betaSink)*ones(1,151);
ntype(1) = -1;

W = zeros(1,numcells);
for n = 1:length(nodes)
    
    % Minimum distance between nodes and atractors/repulsors  
    D = distances(Gnn,nodes(n));
    
    % Energy landscape 
    W = W + ntype(n)*100./(D+1).^0.5; 
    
end

W = (W-max(W));
Pot = exp(W);
Pot = Pot/sum(Pot);

figure('position',[0,0 1800 900])
subplot(1,2,1)
scatter(umap2(:,1),umap2(:,2),[],log(Pot))
colorbar
title('Distribution of log(Differentiation Potential)')
axis off
subplot(1,2,2)
plot(log(sort(Pot)))
title('log(Differentiation Potential)')

filename = strcat({'../../../output/human/figures/log_differentiation_potential.png'});
saveas(gcf, char(filename))

%% Stochastic similations

tic;
s0 = root;
numiter = 2000;
numsimcells = 2000;
states = simcell_potential_1(W,K,s0,numiter,numsimcells);
toc;

%% Plotting simulated cell

mycell = states(:,3);
c2 = colorscatter(mycell,2);

figure()
hold on
scatter(umap2(:,1),umap2(:,2),10,colcluster,'filled')
alpha(.2)
plot(umap2(mycell,1),umap2(mycell,2),'k-');
scatter(umap2(mycell,1),umap2(mycell,2),10,c2,'filled')
axis off 

filename = strcat({'../../../output/human/figures/simulated_cell.png'});
saveas(gcf, char(filename))

%% Average gene expression profiles over stochastic trajectories

numcelltypes = length(clusters);
numend = zeros(numcelltypes,1);

meanG  = cell(numcelltypes,1); % Mean
meanG2 = cell(numcelltypes,1); % Mean square
varG   = cell(numcelltypes,1); % Variance
stdG   = cell(numcelltypes,1); % Standard deviation
semG   = cell(numcelltypes,1); % Standard error of the mean

for i = 1:numcelltypes
    meanG{i}  = zeros(numiter,numgenes); 
    meanG2{i} = zeros(numiter,numgenes); 
end

for c = 1:numsimcells 
    fprintf('%i\n',c);
    ID         = mapcellID(char(cellID(states(numiter,c)))); 
    meanG{ID}  = meanG{ID} + G(states(:,c),:);
    meanG2{ID} = meanG2{ID} + G(states(:,c),:).^2;
    numend(ID) = numend(ID) + 1;
end

for i = 1:numcelltypes
    meanG{i}  = meanG{i}/numend(i);
    meanG2{i} = meanG2{i}/numend(i);
    varG{i}   = meanG2{i}-meanG{i}.^2;
    stdG{i}   = sqrt(varG{i});   
    semG{i}   = stdG{i}/sqrt(numend(i));
end

%% Fate probabilities

pbb = cell(numcelltypes,1);

for i = 1:numcelltypes
    pbb{i} = zeros(numcells,1); 
end
countcell = zeros(numcells,1);
for c = 1:numcells
    fprintf('%i\n',c);
    for t = 1:numsimcells
        ID = mapcellID(char(cellID(states(numiter,t)))); 
        if ismember(c, states(:,t)) == 1
            pbb{ID}(c,1) = pbb{ID}(c,1) + 1;
            countcell(c) = countcell(c) + 1;
        else            
            pbb{ID}(c,1) = pbb{ID}(c,1) + 0;             
        end        
    end    
end

pbb_a = pbb{7}./(countcell+1);
pbb_b = pbb{8}./(countcell+1);
pbb_e = pbb{9}./(countcell+1);

%% Plot distribution of fate probabilities 

figure()
scatter(umap2(:,1),umap2(:,2),5,pbb_a,'filled');
CT=cbrewer('seq', 'YlGnBu', 20);
colormap(CT)
alpha(.6)
set(gca,'FontSize',25)
axis off

filename = strcat({'../../../output/human/figures/alpha_probabilities.png'});
saveas(gcf, char(filename))

figure()
scatter(umap2(:,1),umap2(:,2),5,pbb_b,'filled');
CT=cbrewer('seq', 'YlGn', 20);
colormap(CT)
alpha(.6)
set(gca,'FontSize',25)
axis off

filename = strcat({'../../../output/human/figures/beta_probabilities.png'});
saveas(gcf, char(filename))

figure()
scatter(umap2(:,1),umap2(:,2),5,pbb_e,'filled');
CT=cbrewer('seq', 'YlOrRd', 20);
colormap(CT)
alpha(.4)
set(gca,'FontSize',25)
axis off

filename = strcat({'../../../output/human/figures/ec_probabilities.png'});
saveas(gcf, char(filename))

%% Transcription factor activities CHECKKK!!! Read old file --> check order of TFs 

Act     = importdata('../../../output/human/regularized_activities.csv'); % [factorsxcells] 
A       = Act.data';
TFID    = Act.textdata(:,1);
mapTFID = containers.Map(TFID,1:length(TFID));
numfac  = size(A,2);

%% Average transcription factor activity profiles over stochastic trajectories

numcelltypes = length(clusters);
numend = zeros(numcelltypes,1);

meanTF  = cell(numcelltypes,1); % Mean
meanTF2 = cell(numcelltypes,1); % Mean square
varTF   = cell(numcelltypes,1); % Variance
stdTF   = cell(numcelltypes,1); % Standard deviation
semTF   = cell(numcelltypes,1); % Standard error of the mean

for i = 1:numcelltypes
    meanTF{i}  = zeros(numiter,numfac); 
    meanTF2{i} = zeros(numiter,numfac); 
end

for c = 1:numsimcells 
    fprintf('%i\n',c);
    ID          = mapcellID(char(cellID(states(numiter,c))));
    meanTF{ID}  = meanTF{ID} + A(states(:,c),:);
    meanTF2{ID} = meanTF2{ID} + A(states(:,c),:).^2;
    numend(ID)  = numend(ID) + 1;
end

for i = 1:numcelltypes
    meanTF{i}  = meanTF{i}/numend(i);
    meanTF2{i} = meanTF2{i}/numend(i);
    varTF{i}   = meanTF2{i}-meanTF{i}.^2;
    stdTF{i}   = sqrt(varTF{i}); 
    semTF{i}   = stdTF{i}/sqrt(numend(i));   
end

%% Plotting average gene expression

markers = {'NEUROG3','FEV','PDX1','NKX6-1','ARX','LMX1A'};

%figure('position',[0 0 600 1000])
figure()
for i =1:6
    subplot(2,3,i)
    trajectory_plot_p_h(markers{i},'gene',mapgeneID,meanG,semG,mapTFID,meanTF,semTF)
end

%% Plotting average activity

markers_f = {'ISL1','PDX1','OTP_PHOX2B_LHX1_LMX1A_LHX5_HOXC4','NKX2-2',...
    'NOTO_VSX2_DLX2_DLX6_NKX6-1','PAX6','RFX3_RFX2','ALX1_ARX'};

%figure('position',[0 0 600 1000])
figure()
for i =1:8
    subplot(2,4,i)
    trajectory_plot_p_h(markers_f{i},'tf',mapgeneID,meanG,semG,mapTFID,meanTF,semTF)
end

%% Export data to perform differential motif activity analysis
% This file is read by the script differential_motif_activity_analysis.m

% 7: Alpha, 8: Beta, 9: EC
todel = [1;2;3;4;5;6];

meanG(todel)  = [];
varG(todel)   = [];
semG(todel)   = [];
meanTF(todel) = [];
varTF(todel)  = [];
semTF(todel)  = [];

filename = '../../../output/human/dynamic_trajectories.mat';
save(filename,'states','meanG','varG','semG','mapgeneID',...
    'meanTF','varTF','semTF','mapTFID','TFID')
%%







