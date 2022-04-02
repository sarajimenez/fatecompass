%{
Differential motif activity analysis to identify lineage-specific TFs

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

%% 1. z-score

% Reading data
actDist = load('../../../output/human/bootstrap_distribution_activities.mat');
actDist = actDist.mydist;

% Selecting dimmensions
numfac = length(actDist);
numcells = length(actDist{1});

% Computing z-score
zscore = zeros(numfac,1);
for f = 1:numfac
    
    mu = zeros(numcells,1);
    sigma = zeros(numcells,1);
    for c = 1:numcells
        
        mu(c) = full(mean(actDist{f}(c,:)));
        sigma(c) = full(std(actDist{f}(c,:)));
        
    end
    
    zscore(f) = sqrt(1/numcells * sum((mu./sigma).^2));
    
end

figure()
hold on 
box on
[y,x] = ksdensity(zscore);
plot(x,y,'k','LineWidth',2)
set(gca,'FontSize',25)

filename = strcat({'../../../output/human/figures/zscore_distribution.png'});
saveas(gcf, char(filename))

%% 2. Variability over time 

% Reading data
dynamics = load('../../../output/human/dynamic_trajectories.mat');

meanG  = dynamics.meanG;
varG   = dynamics.varG;
semG   = dynamics.semG;
meanTF = dynamics.meanTF;
varTF  = dynamics.varTF;
semTF  = dynamics.semTF;

% Selecting dimmensions
numcelltypes = length(meanG);

% Standar deviation over time 
stdTF_t = cell(numcelltypes,1);
 
for i = 1:numcelltypes
    stdTF_t{i} = zeros(numfac,1); 
end
 
 
for i = 1:numcelltypes   
    for j = 1:numfac       
        stdTF_t{i}(j,1) = std(meanTF{i}(1:1000,j));        
    end    
end
 
figure()
hold on 
box on
[y,x] = ksdensity(stdTF_t{1});
plot(x,y,'LineWidth',2,'color',A_c)
[y,x] = ksdensity(stdTF_t{2});
plot(x,y,'LineWidth',2,'color',B_c)
[y,x] = ksdensity(stdTF_t{3});
plot(x,y,'LineWidth',2,'color',EC_c)
legend({'Alpha','Beta','EC'})
legend boxoff
hold off

filename = strcat({'../../../output/human/figures/std_over_time_distribution.png'});
saveas(gcf, char(filename))

%% 3. Dynamic cross-correlation between mRNA expression and TF activity

mapgeneID = dynamics.mapgeneID;
TFID      = dynamics.TFID;

CC = cell(numcelltypes,1);
tlag = cell(numcelltypes,1);
 
for i = 1:numcelltypes
    
    for j = 1:numfac
        
        tf = split(TFID{j},"_"); 
        
        if length(tf) == 1
            CC{i}{j,1} = cell(1,1); 
            tlag{i}{j,1} = cell(1,1);
        else
            CC{i}{j,1} = cell(length(tf),1); 
            tlag{i}{j,1} = cell(length(tf),1);            
        end
        
    end
        
end
 
for i = 1:numcelltypes   
    
    for j = 1:numfac
        
        tf = split(TFID{j},"_"); 
        
        if length(tf) == 1 % in this case the motif will always be expressed
        
            [cc,laga] = xcorr((meanTF{i}(:,j)-mean(meanTF{i}(:,j))),...
                (meanG{i}(:,mapgeneID(tf{1}))-mean(meanG{i}(:,mapgeneID(tf{1})))),'coeff'); 
            [maxCC,Ilaga] = max(cc);
            tlag{i}{j,1}(1,1) = num2cell(laga(Ilaga));
            CC{i}{j,1}(1,1) = num2cell(maxCC);
        
        else % in this case some of the tf in the motif may not be expressed 
            
            k = 1; 
            while k <= length(tf)
                
                if isKey(mapgeneID,tf{k}) == 1
                    
                    [cc,laga] = xcorr((meanTF{i}(:,j)-mean(meanTF{i}(:,j))),...
                        (meanG{i}(:,mapgeneID(tf{k}))-mean(meanG{i}(:,mapgeneID(tf{k})))),'coeff'); 
                    [maxCC,Ilaga] = max(cc);
                    tlag{i}{j,1}(k,1) = num2cell(laga(Ilaga));
                    CC{i}{j,1}(k,1) = num2cell(maxCC);
                    
                else 
                    
                    tlag{i}{j,1}(k,1) = num2cell(NaN);
                    CC{i}{j,1}(k,1) = num2cell(NaN);%{'NA'};
                    
                end
                
                k = k + 1;
                
            end
 
        end
        
    end 
    
end

%% Format results of differential motif activity analysis on a table 
  
results = table(TFID,stdTF_t{1}, stdTF_t{2}, stdTF_t{3}, zscore);
 
r = 1;
for j = 1:numfac
    
    tf = split(TFID{j},"_"); 
    
    if length(tf) == 1
        results_all(r,:) = repelem(results(j,:), length(tf), 1);
        tfs(r,1) = tf;
        CCa(r,1) = CC{1}{j,1};
        tlaga(r,1) = tlag{1}{j,1};
        CCb(r,1) = CC{2}{j,1};
        tlagb(r,1) = tlag{2}{j,1};     
        CCe(r,1) = CC{3}{j,1};
        tlage(r,1) = tlag{3}{j,1};
    else
        results_all(r:(r+length(tf)-1),:) = repelem(results(j,:), length(tf), 1);
        tfs(r:(r+length(tf)-1),1) = tf;
        CCa(r:(r+length(tf)-1),1) = CC{1}{j,1};
        tlaga(r:(r+length(tf)-1),1) = tlag{1}{j,1};
        CCb(r:(r+length(tf)-1),1) = CC{2}{j,1};
        tlagb(r:(r+length(tf)-1),1) = tlag{2}{j,1};
        CCe(r:(r+length(tf)-1),1) = CC{3}{j,1};
        tlage(r:(r+length(tf)-1),1) = tlag{3}{j,1};
    end
    r = r + length(tf); 
end
      
varnames = {'motif_family','TF','std_over_alpha_traj','std_over_beta_traj',...
    'std_over_ec_traj','cc_over_alpha_traj','time_max_cc_over_alpha_traj',...
    'cc_over_beta_traj','time_max_cc_over_beta_traj','cc_over_ec_traj',...
    'time_max_cc_over_ec_traj','z_score'};

results_complete = table(table2cell(results_all(:,1)),tfs,...
    table2cell(results_all(:,2)),table2cell(results_all(:,3)),...
    table2cell(results_all(:,4)),...
    CCa,tlaga,CCb,tlagb,CCe,tlage,...
    table2cell(results_all(:,5)),'VariableNames',varnames);
 
writetable(results_complete,...
    '../../../output/human/differential_motif_activity_analysis.txt',...
    'Delimiter', '\t')

%%