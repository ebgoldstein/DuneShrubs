%this is a model of Dunes, Shrubs and TWL started by EBG in 7/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IRBR - 21 Sep 18
% Version includes:
%     + Island width dimension
%     + Seed dispersion in random radial direction
%     + Conversion of shrub age variable to a percent cover variable


% tew dew and open Qs
%-add in VCR TWL data/storms
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
tic
%%%%%%%%%%%%%%%%%%
%TIME
TMAX=30 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the model computational domain, in meters
DomainLength=500;
DomainWidth=50;
DuneDomain=zeros(TMAX,DomainLength);
ShrubDomainAll=zeros(TMAX,DomainLength,DomainWidth);
ShrubDomainFemale=ShrubDomainAll;
ShrubDomainMale=ShrubDomainAll;
ShrubPercentCover=ShrubDomainAll;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the Dune parameters:
%initial elevation of the dune domain
DuneDomain(1,:)=ones(1,DomainLength)*0.1;

%Maximum dune height
Dmax=5;

%Percent cover change
PC = [4 8 10 15 15 20 35 80 100 100];

%Dune intrinsic growth rate
%houser et al 0.05 ? r ? 0.55
rmin=0.05;
rmax=0.55;
growthrate= rmax + (rmax-rmin).*rand(1,DomainLength);

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
ShrubDomainFemale(1,round(DomainLength/2),round(DomainWidth/2))=1;
% ShrubDomainFemale(1,round(DomainLength*0.33),1)=1;
% ShrubDomainFemale(1,round(DomainLength*0.66),1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yearly time steps
fprintf('\nTime step: ');

for t=2:1:TMAX
    
    % Print time step to screen
    if t > 1
        for j = 0:log10(t-1)
            fprintf('\b');
        end
    end
    fprintf('%d', t);
    
    %Dune Growth
    G=growthrate.*(DuneDomain(t-1,:)).*(1-(DuneDomain(t-1,:)/Dmax));
    DuneDomain(t,:)=G+DuneDomain(t-1,:);
    
    %Shrubs
    %age the shrubs in years
    for k = 1:DomainWidth %Loop through each row of island width (i.e. from ocean to mainland side of island)
        growF=(ShrubDomainFemale(t-1,:,k)>0) & (ShrubDomainFemale(t-1,:,k)<10);
        growM=(ShrubDomainMale(t-1,:,k)>0) & (ShrubDomainMale(t-1,:,k)<10);
        ShrubDomainFemale(t,:,k)=ShrubDomainFemale(t-1,:,k)+ growF;
        ShrubDomainMale(t,:,k)=ShrubDomainMale(t-1,:,k)+ growM;
        ShrubDomainAll(t,:,k)=ShrubDomainMale(t,:,k)+ShrubDomainFemale(t,:,k);

        % Find percent cover for each cell
        for j = 1:(length(ShrubDomainAll)) %Loop through each shrub cell in domain
            a = ShrubDomainAll(t,j,k); %Get age of shrub in cell
            if a > 10
                a = 10;
            end
            
            if a == 0
                ShrubPercentCover(t,j,k) = 0; %Zero percent cover if age is zero
            else
                p = PC(a); %Get corresponding percent cover for its age
                ShrubPercentCover(t,j,k) = p;
            end
        end
    end
    
    for k = 1:DomainWidth %Loop through each row of island width (i.e. from ocean to mainland side of island)
        % for all cells with a shrub
        FemaleShrubs = squeeze(ShrubDomainFemale(:,:,k));
        I = find(FemaleShrubs(t,:)>=TimeFruit);
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
                        %   -there is no shrub at the receiving cell (male or
                        %   female)
                        %   -the receiving cell has a tall enough fronting dune
                        if  targetX > 0 && targetX <= DomainLength && targetY > 0 && targetY <= DomainWidth %If target point is within island domain
                            if ShrubDomainAll(t,targetX,targetY)==0 && DuneDomain(t,targetX)>=Dshrub % If no shrubs exist there and dune is tall enough
                                %decide if the tree wll be a female or male
                                if rand > Female
                                    ShrubDomainFemale(t,targetX,targetY)=1;
                                else
                                    ShrubDomainMale(t,targetX,targetY)=1;
                                end
                            end
                        end
                    end
                 end
            end
        end   
    %Storms
    end
end

ShrubDomainAll=ShrubDomainMale+ShrubDomainFemale;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(3,1,1), imagesc(DuneDomain);figure(gcf);
title('Dune Height');
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')

ShrubRow1 = squeeze(ShrubDomainAll(:,:,1)); % Plots only the front row of the island over time
subplot(3,1,2), imagesc(ShrubRow1);figure(gcf);
title('Shrub Age');
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')

PCRow1 = squeeze(ShrubPercentCover(:,:,1)); % Plots only the front row of the island over time
subplot(3,1,3), imagesc(PCRow1);figure(gcf);
title('Shrub Percent Cover');
colorbar
axis xy
ylabel('t')
xlabel('alongshore distance')

% Create movie frames of each time-step
close all
DuneShrubMovie(ShrubDomainAll,TMAX,1);
close all
DuneShrubMovie(ShrubPercentCover,TMAX,2);


% figure
% plot(CompDomain(5,:))

%test as a function of dispersal kernal function...
%can be matched to the linear invadsion rate of shrub front on Hog:
%4km/ 26 years ? 153m/yr calculated from Fig1 of Zinnert 2011)

fprintf('\nElapsed time: %.2f\n', toc); % Print elapsed time of program
