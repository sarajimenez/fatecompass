%{
Simulation of stochastic trajectories during pancreas development
Data: Bastidas-Ponce et. al 2019 
System: Mouse in-vivo at E15.5

Author: Sara Jimenez
Email: sarajcorrea@gmail.com
Molina's lab and Gradwohl's lab
IGBMC | Strasbourg University
%}


clc; clear; close all; 

% CAREFUL: this must be an absolute path --> change accordingly
fcnpath = '/Users/jimenezs/Sara/fatecompass/reproducibility/functions';
userpath(fcnpath)

%% Import pre-processed data

datapath = '../../data/mouse/';

% Reading data
umap2  = importdata(strcat(datapath,'umap_2.txt'));
velo2  = importdata(strcat(datapath,'velo_umap_2.txt'));
umap10 = importdata(strcat(datapath,'umap_10.txt'));
velo10 = importdata(strcat(datapath,'velo_umap_10.txt'));
cellID = importdata(strcat(datapath,'cell_id.txt'));
GE     = importdata(strcat(datapath,'normalized_gene_expression.tsv'));
G      = GE.data'; % [cellsxgenes]

% Selecting dimmensions
[numcells,numgenes] = size(G);

% Containers to map cell-types and gene names
clusters = {'Ductal','Ngn3 low EP','Ngn3 high EP', 'Pre-endocrine', ...
    'Beta', 'Alpha','Delta', 'Epsilon'};
mapcellID = containers.Map(clusters,1:length(clusters));
mapgeneID = containers.Map(GE.textdata,1:length(GE.textdata)); 

%% Figure parameters

Acol = [0 121 185]/255;
Bcol = [165 225 127]/255;

colcluster = strcmp(cellID,'Ductal')*[153 191 153]/255 ...
    + strcmp(cellID,'Ngn3 low EP')*[255 153 181]/255 ...
    + strcmp(cellID,'Ngn3 high EP')*[255 189 95]/255 ...
    + strcmp(cellID,'Pre-endocrine')*[255 118 0]/255 ...
    + strcmp(cellID,'Beta')*[165 225 127]/255 ...
    + strcmp(cellID,'Alpha')*[0 121 185]/255 ...
    + strcmp(cellID,'Delta')*[114 56 159]/255 ...
    + strcmp(cellID,'Epsilon')*[207 176 217]/255;

%% Nearest neighbor graph representing phenotypic manifold

numnb = 10;

Dd = pdist2(umap10, umap10);
Sv = sort(Dd);
NN = Sv(1:numnb,:);

%% Modeling transition probabilities using a Markov process

% Find optimal parameters 
% Sets dt such that on average the number of nearest neighbors can be reached
% Sets D such that the average number of connection is 2*numnb 

cutoff = 0.01;
D0 = 0.05; % initial condition 

dt =  mean(NN(:))/mean(vecnorm(velo10')); 
Dv = pdist2(umap10, umap10+velo10*dt); 
avercon = @(D) average_number_connections(Dv,D,cutoff)-2*numnb;
options = optimset('Display','iter');
D = fzero(avercon,D0,options);
fprintf('D = %f dt = %f\n',D,dt);

% Calculating normalized transition probabilities
Pv = exp(-(Dv.^2)/D);
Pv = Pv./sum(Pv);
Pv = Pv./sum(Pv);
% Propagator
P = Pv;

% Transitions in cells for simulation
W = cell(numcells,1);
N = cell(numcells,1);
for i = 1:numcells
    N{i} = find(P(:,i)); % indixes of all the posible transitions for cell i
    W{i} = P(N{i},i); % weigth of all the posible transitions for cell i 
    W{i} = W{i}/sum(W{i}); % normalized weigth of all the posible transitions for cell i 
end

%% Plots to check if D and dt make sense 

figure('position',[0 0 1200 500])

subplot(2,4,1)
hold on
box on
[x, y] = ksdensity(NN(:));
plot(y, x, 'LineWidth', 2)
[x, y] = ksdensity(vecnorm(velo10')*dt);
plot(y, x, 'LineWidth', 2)
legend({'Distance to nearest nb','Distance traveled dt*|V|'})
set(gca,'FontSize',15)
legend boxoff

subplot(2,4,5)
ksdensity(sum(P > cutoff))
legend({'Degree distribution'})
set(gca,'FontSize',15)
legend boxoff

subplot(2,4,[2 3 4 6 7 8])
hold on
quiver(umap2(:,1),umap2(:,2),velo2(:,1),velo2(:,2),'color',[145 145 145]/255)
scatter(umap2(:,1),umap2(:,2),20,colcluster,'filled')
alpha(.2)

%Alpha
quiver(umap2(P(:,7)>cutoff,1),umap2(P(:,7)>cutoff,2),velo2(P(:,7)>cutoff,1),velo2(P(:,7)>cutoff,2),'k','LineWidth',2)
scatter(umap2(P(:,7)>cutoff,1),umap2(P(:,7)>cutoff,2),20,'r','filled')
quiver(umap2(7,1),umap2(7,1),velo2(7,1),velo2(7,1),'k','LineWidth',2)
scatter(umap2(7,1),umap2(7,2),80,[0 121 185]/255,'filled')

%Beta
quiver(umap2(P(:,17)>cutoff,1),umap2(P(:,17)>cutoff,2),velo2(P(:,17)>cutoff,1),velo2(P(:,17)>cutoff,2),'k','LineWidth',2)
scatter(umap2(P(:,17)>cutoff,1),umap2(P(:,17)>cutoff,2),20,'r','filled')
quiver(umap2(17,1),umap2(17,1),velo2(17,1),velo2(17,1),'k','LineWidth',2)
scatter(umap2(17,1),umap2(17,2),80,[165 225 127]/255,'filled')

% Ngn3 high
quiver(umap2(P(:,3592)>cutoff,1),umap2(P(:,3592)>cutoff,2),velo2(P(:,3592)>cutoff,1),velo2(P(:,3592)>cutoff,2),'k','LineWidth',2)
scatter(umap2(P(:,3592)>cutoff,1),umap2(P(:,3592)>cutoff,2),20,'r','filled')
quiver(umap2(3592,1),umap2(3592,1),velo2(3592,1),velo2(3592,1),'k','LineWidth',2)
scatter(umap2(3592,1),umap2(3592,2),80,[255 189 95]/255,'filled')
axis off

filename = strcat({'../../../output/mouse/figures/fit_D_dt.png'});
saveas(gcf, char(filename))

%% Stochastic similations 

tic;
s0 = 336; % Sox9+ bipotent progenitor
%s0 = 3592; % Ngn3 high
numiter = 2000; 
numsimcells = 5000;
states = simcell(W,N,s0,numiter,numsimcells);
toc;

%% Plot simulated cell

mycell = states(:,3);
c2 = colorscatter(mycell,2);

figure()
hold on 
quiver(umap2(:,1),umap2(:,2),velo2(:,1),velo2(:,2),'color',[145 145 145]/255)
scatter(umap2(:,1),umap2(:,2),20,colcluster,'filled')
alpha(.2)
plot(umap2(mycell,1),umap2(mycell,2),'k-');
scatter(umap2(mycell,1),umap2(mycell,2),10,c2,'filled')
axis off

filename = strcat({'../../../output/mouse/figures/simulated_cell.png'});
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

pbb_a = pbb{6}./(countcell+1);
pbb_b = pbb{5}./(countcell+1);

%% Plot distribution of fate probabilities

figure()
scatter(umap2(:,1),umap2(:,2),20,pbb_a,'filled');
CT=cbrewer('seq', 'YlGnBu', 20);
colormap(CT)
set(gca,'FontSize',25)
axis off

filename = strcat({'../../../output/mouse/figures/alpha_probabilities.png'});
saveas(gcf, char(filename))

figure()
scatter(umap2(:,1),umap2(:,2),20,pbb_b,'filled');
CT=cbrewer('seq', 'YlGn', 20);
colormap(CT)
set(gca,'FontSize',25)
axis off

filename = strcat({'../../../output/mouse/figures/beta_probabilities.png'});
saveas(gcf, char(filename))

%% Transcription factor activities 

Act     = importdata('../../../output/mouse/regularized_activities.csv'); % [factorsxcells] 
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

markers = {'Neurog3','Fev','Pdx1','Nkx6-1','Arx','Irx2'};

figure('position',[0 0 600 1000])

for i =1:6
    subplot(3,2,i)
    trajectory_plot(markers{i},'gene',mapgeneID,meanG,semG,mapTFID,meanTF,semTF)
end

%% Plotting average activity 

markersTF = {'Neurod1','Nkx2-2','Hoxb8_Pdx1',...
    'Nkx6-1_Evx1_Hesx1','Hoxc4_Arx_Otp_Esx1_Phox2b','Irx2_Irx3','Pax4','Arntl_Tfe3_Mlx_Mitf_Mlxipl_Tfec'};

figure('position',[0 0 600 1000])

for i =1:8
    subplot(4,2,i)
    trajectory_plot(markersTF{i},'tf',mapgeneID,meanG,semG,mapTFID,meanTF,semTF)
end

%% Export data to perform differential motif activity analysis
% This file is read by the script differential_motif_activity_analysis.m

% 6: Alpha, 5: Beta
todel = [1;2;3;4;7;8];

meanG(todel)  = [];
varG(todel)   = [];
semG(todel)   = [];
meanTF(todel) = [];
varTF(todel)  = [];
semTF(todel)  = [];

filename = '../../../output/mouse/dynamic_trajectories.mat';
save(filename,'states','meanG','varG','semG','mapgeneID',...
    'meanTF','varTF','semTF','mapTFID','TFID')
%%

