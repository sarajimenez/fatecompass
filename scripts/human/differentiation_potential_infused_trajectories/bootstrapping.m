% v4: Smooth regularization of activities --> correction powering ONLY M
% and changing neighborhood graph added line 97 and changed line 101

%{
Smooth regularization - Linear model to estimate TF activities
Data: Bastidas-Ponce et. al 2019 
System: Mouse pancreata - endocrine cell at E15.5

Sara Jimenez
Gradwohl's lab and Molina's lab 
IGBMC

%}

clc; clear; close all; 
newpath = '/Users/jimenezs/Sara/CellFate_Diffusion/functions';
userpath(newpath)

%% 
tic;

U2 = importdata('../../Data/potential/umap2.txt');
U10 = importdata('../../Data/potential/umap10.txt'); % 30 was the number of neighboors used to compute the umap

E = importdata('../../Data/potential/GE_new.tsv');
E = E.data;
E = E'; % cells x genes
N = importdata('../../Data/potential/BSperFamily.csv');
N = N.data; % genes x factors
toc;

% Getting dimensions of the problem
tic;
E = sparse(E);
N = sparse(N);
[numcells,numgenes] = size(E);
[~,numfac]   = size(N);

%% Normalization for activity estimation
tic;

E = E - (1/numgenes)*sum(E,2);
E = E - (1/numcells)*sum(E,1);

Nm = (1/numgenes)*sum(N);
N = N - Nm;

fprintf('-- Expression matrix cell- and gene- normalized --\n');
fprintf('-- Binding site matrix motif-normalized  --\n');
toc;

%% Creating the network (This is not very efficient. You can change it with the way you do it)

tic
numberconnections = 51;
I = knnsearch(U10,U10,'K',numberconnections);
M = zeros(numcells);
for i = 1:numcells
    M(i,I(i,2:end)) = 1;%1./pdist2(E(i,:),E(I(i,2:end,:),:)).^2; % The inverse distance between cell i and its neighbours 
end

M = M + 10*eye(size(M)); % Self transitions are important
M = M./sum(M,2);
M = M';

fprintf('-- Network computed and row-normalized --\n');
%plot(digraph(M),'XData',U2(:,1),'YData',U2(:,2))

M = sparse(M);

toc;

%%
tic

numsamples  = 5;
Est_dist    = cell(numsamples,1);  

for i = 1:numsamples
        
    genestr = rand(numgenes,1)<0.8;
    Ntr = N(genestr,:);
    Etr = E(:,genestr);

    Atr = lsqminnorm(Ntr,Etr'); % factors x cells 

    Est_dist{i} = M^4 * Atr'; %cells x factors
    
    fprintf('%i\n',i)
        
end

toc

%% Formating the distributions

tic

mydist = cell(numfac,1);

for f = 1:numfac
    
    for i = 1:numsamples
        
        mydist{f}(:,i) = Est_dist{i}(:,f);
        
    end
      
end

toc

%% Load distribution of the estimate with 1000 samples

mydist = load('/Volumes/gradwohl/USERS/User_Sara/Published_scRNA/Veres_2019/Analysis_S5_timecourse/activities/mydist.mat');
mydist = mydist.mydist;

%%

myfac = 16;
mycell = 2000;

figure()
subplot(1,2,1)
histogram(mydist{myfac}(mycell,:),100);
title('Distribution of the sample mean')

subplot(1,2,2)
ksdensity(mydist{myfac}(mycell,:));
title('Distribution estimate of the sample mean')



%% z-score
tic;

zscore = zeros(numfac,1);
for f = 1:numfac
    
    mu = zeros(numcells,1);
    sigma = zeros(numcells,1);
    for c = 1:numcells
        
        mu(c) = full(mean(mydist{f}(c,:)));
        sigma(c) = full(std(mydist{f}(c,:)));
        
    end
    
    zscore(f) = sqrt(1/numcells * sum((mu./sigma).^2));
    
end

toc;

%% save in the regularization folder
% 
dlmwrite('output/z_score_100.tsv', zscore, 'delimiter','\t')

%dlmwrite('output/Areg.csv', Areg, 'delimiter',',')



















