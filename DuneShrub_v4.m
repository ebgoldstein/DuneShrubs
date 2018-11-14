% A model of Dunes, Shrubs and TWL 

% started by EBG in 7/2017; IRBR 10-11/18
% Version 4 includes:
%     + NO elevation
%     + Storms at frequency and magnitude according to Hog Island TWL time series
%     + Dune erosion following Long et al. (2014) data and maximum slope angle enforcement
%     + Simplistic shrub mortality

% OTHER Qs? 
% - do shrubs die with storm? (probably not if mature, i.e., VCU papers..)
% - fecundity as a fn of maturity (in this case, percent cover)?
% - Seed Survival; i.e., predation by mice and rats?
% - longevity of trees in cell before they stop producing seeds?
% - seed crop interval?
% - germination rate?
% - mortality?

clear all
close all
clc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET PARAMETERS       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME
TMAX = 15;
StormStart = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTATIONAL DOMAIN (meters)
DomainLength = 1000;
DomainWidth = 100;
DuneDomain = zeros(TMAX,DomainLength);
ShrubDomainAll = zeros(TMAX,DomainLength,DomainWidth);
ShrubDomainFemale = ShrubDomainAll;
ShrubDomainMale = ShrubDomainAll;
ShrubPercentCover = ShrubDomainAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DUNE
% Initial elevation of the dune domain
% DuneDomain(1,:)=ones(1,DomainLength)*0.1;
DuneDomain(1,:) = ones(1,DomainLength) .* (1.1 + (-0.1 + (0.1 - (-0.1))*rand(1,DomainLength)));

% Maximum dune height
Dmax = 3; % Originally 5

% Dune intrinsic growth rate
% Houser et al 0.05 ? r ? 0.55
rmin = 0.05;
rmax = 0.55;
growthrate = rmax + (rmax-rmin).*rand(1,DomainLength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STORM

% Number Per Year
A_storm = 52.398; % For Weibull (temp)
B_storm = 4.5301;   

% Total Water Level
mu_storm = 1.1075; %1.1075
sigma_storm = 0.3353; %0.3353

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SHRUB

% Shrub Growth from JCZ
% Minimum dune height for shrub growth
Dshrub = 2;

% Percent cover change (year 1-10)
PC = [4 8 10 15 15 20 35 80 100 100];

% Seeds produced per cell per year (Fecundity)
% Seedmin=1000; %???
% Seedmax=10000; % from Kwit et al 2004 Oecologia
Seedmin = 100; % Temp
Seedmax = 1000; % Temp

% Seed dispersal distance
% Dows 2016 has the data
% i took the 0 thicket data from 0-2m and then the grassland data from 50-300m. then fit a lognormal distribution (better than exponential) to data.
% need to use fn 'lognrnd' to get dispersal distances. i.e., DispDist = lognrnd(mu,sigma)
mu = -0.721891;
sigma = 1.5;

% Germination rate; 60% from Young et al 1994
GermRate = 0.6;

% Survival rate (from those reported in Lazarus and McGill)
%SurvRatemin=10e-5;
%SurvRatemax=0.02;
SurvRatemin = GermRate;
SurvRatemax = GermRate;

% Trees take 4-6 years to fruit; email from Julie
TimeFruit = 5;

% Trees are 50/50 male/female; Hokkanen thesis 
Female = 0.5;

% Initiate Shrubs
ShrubDomainFemale(1,round(DomainLength/2),round(DomainWidth/2)) = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stormCount = zeros(1,2); % Stores year and TWL of every storm 

fprintf('\nTime step: ');

for t = 2:1:TMAX % Yearly time steps
    
    % Print time step to screen
    if t > 1
        for j = 0:log10(t-1)
            fprintf('\b');
        end
    end
    fprintf('%d', t);

    %%% Dune Growth
    G = growthrate.*(DuneDomain(t-1,:)).*(1-(DuneDomain(t-1,:)/Dmax));
    DuneDomain(t,:) = G + DuneDomain(t-1,:);
    
    %%% Shrubs
    % Age the shrubs in years
    for k = 1:DomainWidth %Loop through each row of island width (i.e. from ocean to mainland side of island)
        growF = (ShrubDomainFemale(t-1,:,k) > 0) & (ShrubDomainFemale(t-1,:,k) < 10);
        growM = (ShrubDomainMale(t-1,:,k) > 0) & (ShrubDomainMale(t-1,:,k) < 10);
        ShrubDomainFemale(t,:,k) = ShrubDomainFemale(t-1,:,k) + growF;
        ShrubDomainMale(t,:,k) = ShrubDomainMale(t-1,:,k) + growM;
        ShrubDomainAll(t,:,k) = ShrubDomainMale(t,:,k) + ShrubDomainFemale(t,:,k);

        % Find percent cover for each cell
        for j = 1:(length(ShrubDomainAll)) % Loop through each shrub cell in domain
            a = ShrubDomainAll(t,j,k); % Get age of shrub in cell
            if a > 10
                a = 10;
            end
            if a == 0
                ShrubPercentCover(t,j,k) = 0; % Zero percent cover if age is zero
            else
                p = PC(a); % Get corresponding percent cover for its age
                ShrubPercentCover(t,j,k) = p;
            end
        end
    end
    
    %%% Disperse seeds
    for k = 1:DomainWidth % Loop through each row of island width (i.e. from ocean to mainland side of island)
        % For all cells with a shrub
        FemaleShrubs = squeeze(ShrubDomainFemale(:,:,k));
        I = find(FemaleShrubs(t,:)>=TimeFruit);
        numShrubCells=length(I);
        % Determine how many seeds in each cell
        Seedspercell= randi([Seedmin,Seedmax],1,numShrubCells);
        % For each shrub cell, determine survival rate for the cell in this year
        SeedSurvYear= SurvRatemin + (SurvRatemax-SurvRatemin).*rand(numShrubCells,1);

        for i=1:numShrubCells
            % For each shrub cell producing seeds, generate a random # of random numbers each representing a single seed
            randseeds=rand(1,Seedspercell(i));
            % Find how many seeds produced in each cell actually survive
            Survivors=length(randseeds(randseeds<SeedSurvYear(i)));

            % Determine distance, rounding to nearest integer
            if Survivors>0
                DispDist = round(lognrnd(mu,sigma,1,Survivors));
            else
                DispDist=[];
            end

            % If there are seeds to disperse
            if isempty(DispDist)==0
                for j=1:length(DispDist) % Loop through each individual seed to disperse

                    % Create a meshgrid to find coordinates of all points that are dropdistance from origin 
                    if DispDist(j)>0
                        [X,Y] = meshgrid(1:(DispDist(j)*2+1),1:(DispDist(j)*2+1));
                        originX = I(i); % Sets coordinates of plant where seed is coming from
                        originY = k;
                        matOriginX = round(length(X)/2);
                        matOriginY = matOriginX;
                        distMat = round(sqrt((X-matOriginX).^2 + (Y-matOriginY).^2)); % Find the distance from origin to every other point on island
                        [row,col] = find(distMat == DispDist(j)); % Select coordinates of each point on island that is dropdistance away from origin
                        
                        % Randomly selct one of those points as the target point - this means equal probability of dispersal in every direction
                        if isempty(col)
                            targetX = originX;
                            targetY = originY;
                        else
                            xx = randi(length(col));
                            matTargetX = col(xx);
                            matTargetY = row(xx);
                            
                            targetX = originX + (matTargetX-matOriginX);
                            targetY = originY + (matTargetY-matOriginY);
                        end
                        
                        % Put a plant in the ground if 
                        %   -the dropdistance>0, 
                        %   -the target drop location is within the island domain, 
                        %   -there is no shrub at the receiving cell (male or female)
                        %   -the receiving cell has a tall enough fronting dune
                        if  targetX > 0 && targetX <= DomainLength && targetY > 0 && targetY <= DomainWidth ... 
                            && ShrubDomainAll(t,targetX,targetY) == 0 && DuneDomain(t,targetX) >= Dshrub
                            % Decide if the tree wll be a female or male
                            if rand > Female
                                ShrubDomainFemale(t,targetX,targetY) = 1;
                            else
                                ShrubDomainMale(t,targetX,targetY) = 1;
                            end
                        end
                    end
                 end
            end
        end   
    end
    
    %%% Storms
    if t >= StormStart
        % Select number of storms for this time step from weibull distribution
        numstorm = round(wblrnd(A_storm,B_storm));

        if numstorm > 0
            % Select TWL of each storm from normal distribution
            TWL = normrnd(mu_storm, sigma_storm, 1, numstorm);
            for n = 1:length(TWL)
                stormCount = [stormCount; t TWL(n)]; % Store all storm records in array
            end  

            % Overwash/Erosion
            % Reduce dune height for overwashed dunes & kill shrubs
            for n = 1:numstorm
                Low = unifrnd(0.3,0.8); % Long et al. (2014) - 80% of observed change fall within 30-80% of pre-storm dune height
                Dow = find(DuneDomain(t,:) < TWL(n)); % All overwashed dune cells
                OWdist = randi([30 80], length(Dow), 1); % Random overwash penetration distance
                %                 DuneDomain(DuneDomain(t,:) < TWL(n)) = DuneDomain(DuneDomain(t,:) < TWL(n)) * Low; % Lower dunes -- Doesn't work
                
                
                %Remove the relevant slice of the array to opertate on.
                ShrubMortF=squeeze(ShrubDomainFemale(t,:,:));
                ShrubMortM=squeeze(ShrubDomainFemale(t,:,:));
                for d = 1:length(Dow)
                    DuneDomain(t,Dow(d)) = DuneDomain(t,Dow(d)) * Low; % Less efficient - but works
                    % Kill 1-year-old shrubs behind overwashed dunes
                    ShrubMortF(ShrubMortF(Dow(d),1:OWdist(d)) == 1) = 0; % Doesn't seem to work
                    ShrubMortM(ShrubMortM(Dow(d),1:OWdist(d)) == 1) = 0; % Doesn't seem to work
                    % Kill all shrubs behind dunes lowered below threshold height
                end
                %Put the relevant slice of the array back in.
                ShrubDomainFemale(t,:,:)=ShrubMortF;
                ShrubDomainMale(t,:,:)=ShrubMortM;
                
                % Kill all shrubs behind dunes lowered below threshold height
                for d = 1:length(Dow)
                    if DuneDomain(Dow(d)) < Dshrub
                        ShrubDomainFemale(t,Dow(d),1:OWdist(d)) = 0;
                        ShrubDomainMale(t,Dow(d),1:OWdist(d)) = 0;
                    end
                end
                
                % Enforce maximum dune slope angle (37% from Rastetter model, i.e. 0.8 m in height diff betwen 1 m long dune cells)
                for d = 2:(DomainLength) % Loop L to R
                    Lslope = DuneDomain(t,d) - DuneDomain(t,d-1); % Calculate height difference
                    % Subtract sand from taller dune cell and add to adjacent shorter one (assumes sand within dunes is conserved)
                    if Lslope > 0.8
                        sub = (Lslope - 0.8)/2; % Height difference in excess of maximum, divided by two
                        DuneDomain(t,d) = DuneDomain(t,d) - sub; % Lower dune that's too tall
                        DuneDomain(t,d-1) = DuneDomain(t,d-1) + sub; % Raise dune that's too short
                    end
                end
                for d = (DomainLength-1):-1:1 % Loop R to L
                    Rslope = DuneDomain(t,d) - DuneDomain(t,d+1);
                    if Rslope > 0.8
                        sub = (Rslope - 0.8)/2;
                        DuneDomain(t,d) = DuneDomain(t,d) - sub;
                        DuneDomain(t,d+1) = DuneDomain(t,d+1) + sub;
                    end
                end
            end
        end
    end
end    

ShrubDomainAll = ShrubDomainMale + ShrubDomainFemale;
stormCount = stormCount(2:end,:);

% End of model run
fprintf('\nElapsed time: %.2f\n', toc); % Print elapsed time of simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot dunes over time
figure
imagesc(DuneDomain);figure(gcf);
title('Dune Height');
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')



% Create movie frames of each time-step
figure()
DuneShrubMovie(ShrubDomainAll,TMAX,1);
% figure()
% DuneShrubMovie(ShrubPercentCover,TMAX,2);


% Plot storms
% figure()
% scatter(stormCount(:,1),stormCount(:,2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0])



%test as a function of dispersal kernal function...
%can be matched to the linear invadsion rate of shrub front on Hog:
%4km/ 26 years ? 153m/yr calculated from Fig1 of Zinnert 2011)
