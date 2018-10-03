%this is a model of Dunes, Shrubs and TWL started by EBG in 7/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tew dew and open Qs
%-add in VCR TWL data
%-shrub growth (percent cover)
%-expansion rate should be faster, so BI needs to have width...
%OTHER Qs?
%- do shrubs die with storm? (probably not if mature, i.e., VCU papers..)
%-fecundity as a fn of maturity (in this case, percent cover)?
%-Seed Survival; i.e., predation by mice and rats?
%-longevity of trees in cell before they stop producing seeds?
%-seed crop interval?
%-germination rate?
%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%
%TIME
TMAX=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the model computational domain, in meters
DomainSize=1000;
DuneDomain=zeros(TMAX,DomainSize);
ShrubDomainMale=DuneDomain;
ShrubDomainFemale=DuneDomain;
ShrubDomainAll=DuneDomain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the Dune parameters:
%initial elevation of the dune domain
DuneDomain(1,:)=ones(1,DomainSize)*0.1;

%Maximum dune height
Dmax=5;

%Dune intrinsic growth rate
%houser et al 0.05 ? r ? 0.55
rmin=0.05;
rmax=0.55;
growthrate= rmax + (rmax-rmin).*rand(1,DomainSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shrub parameters

%Shrub Growth from Julie's email

%Minimum dune height for shrub growth
Dshrub=2;

%Seeds produced per cell per year (Fecundity), 
Seedmin=1000; %???
Seedmax=10000; %from Kwit et al 2004 Oecologia

%Seed dispersal distance
% Dows 2016 has the data
% i took the 0 thicket data from 0-2m and then the grassland data from
% 50-300m. then fit a lognormal distribution (better than exponential) to
% data.
%  Lognormal distribution
%       mu = -0.721891   [-0.773302, -0.67048]
%    sigma =   1.08424   [1.04909, 1.12184]
% %%%%%%%%%%%%%%%%%%
% load Dows
% load Dowstwo
% pd=fitdist(Dowstwo','Lognormal')
% x_values = 0.1:0.1:500;
% y = pdf(pd,x_values);
% semilogx(x_values,y,'LineWidth',2)
% hold on
% plot(Dows(:,1),Dows(:,2)/sum(Dows(:,2)),'or')
% %%%%%%%%%%%%%%%%%%
% Nathan and Muller-Landau 2000 discuss errors associated with kernal data like this:
% - Assuming source is closest will underestimate dispersal
% - Assuming source is equiprobable for all possible dispersers will overestiamte dispersal distsnce. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%need to use fn 'lognrnd' to get dispersal distances. i.e., DispDist = lognrnd(mu,sigma)
mu = -0.721891;
sigma =   1.5;

%Germination rate; 60% from Young et al 1994;
GermRate=0.6;

%the Survival rate come from those reported in Lazarus and McGill;
% SurvRatemin=10e-5;
% SurvRatemax=0.02;

SurvRatemin=GermRate;
SurvRatemax=GermRate;

%Trees take 4-6 years to fruit; email from Julie
TimeFruit=5;

%Trees are 50/50 male/female; Hokkanen thesis 
Female=0.5;

%OTHER processes???

%Initiate Shrubs
ShrubDomainFemale(1,1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yearly time steps

for t=2:1:TMAX
    
    %Dune Growth
    G=growthrate.*(DuneDomain(t-1,:)).*(1-(DuneDomain(t-1,:)/Dmax));
    DuneDomain(t,:)=G+DuneDomain(t-1,:);
    
    %Shrubs
    
    %age the shrubs in years, doesn't grow them yet in percent cover
    growF=(ShrubDomainFemale(t-1,:)>0) & (ShrubDomainFemale(t-1,:)<10);
    growM=(ShrubDomainMale(t-1,:)>0) & (ShrubDomainMale(t-1,:)<10);
    ShrubDomainFemale(t,:)=ShrubDomainFemale(t-1,:)+ growF;
    ShrubDomainMale(t,:)=ShrubDomainMale(t-1,:)+ growM;
    ShrubDomainAll(t,:)=ShrubDomainMale(t,:)+ShrubDomainFemale(t,:);
    
    % for all cells with a shrub
    I = find(ShrubDomainFemale(t,:)>=TimeFruit);
    numShrubCells=length(I);
    % determine how many seeds in each cell
    Seedspercell= randi([Seedmin,Seedmax],1,numShrubCells);
    % for each shrub cell, determine survival rate for the cell in this
    % year
    SeedSurvYear= SurvRatemin + (SurvRatemax-SurvRatemin).*rand(numShrubCells,1);
    
    for i=1:numShrubCells
        %for each shrub cell producing seeds, generate a random # of
        %random numbers each representing a single seed
        %seeds
        randseeds=rand(1,Seedspercell(i));
        %Find how many seeds produced in each cell actually survive
        Survivors=length(randseeds(randseeds<SeedSurvYear(i)));
        
        % pick a random direction (only 2 directions matter)
        %in a sense the dispersal kernal is one sided, so for now all the
        %seeds will disperse in the along-aisland direction
        
        % determine distance, rounding to nearest integer
        if Survivors>0
            DispDist = round(lognrnd(mu,sigma,1,Survivors));
        else
            DispDist=[];
        end
        
        % if there are seeds to disperse
        if isempty(DispDist)==0
            for j=1:length(DispDist)
                % Put a plant in the ground if 
                %   -the dropdistance>0, 
                %   -the dropdistance less than the length of the island, 
                %   -the there is no shrub at the receiving cell (male or
                %   female)
                %   -the receiving cell has a tall enough fronting dune 
                if DispDist(j)>0 && DispDist(j)< DomainSize-I(i)...
                        && ShrubDomainFemale(t,I(i)+DispDist(j))==0 ...
                        && ShrubDomainMale(t,I(i)+DispDist(j))==0 ...
                        && DuneDomain(t,I(i)+DispDist(j))>=Dshrub
                    %decide if the tree wll be a female or male
                    if rand > Female
                        ShrubDomainFemale(t,I(i)+DispDist(j))=1;
                    else
                        ShrubDomainMale(t,I(i)+DispDist(j))=1;
                    end
                end
            end
        end
    end   
    %Storms
end

ShrubDomainAll=ShrubDomainMale+ShrubDomainFemale;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1), imagesc(DuneDomain);figure(gcf);
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')

subplot(2,1,2), imagesc(ShrubDomainAll);figure(gcf);
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')

% figure
% plot(CompDomain(5,:))

%test as a function of dispersal kernal function...
%can be matched to the linear invadsion rate of shrub front on Hog:
%4km/ 26 years ? 153m/yr calculated from Fig1 of Zinnert 2011)

