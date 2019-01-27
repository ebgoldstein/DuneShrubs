
# A model of Dunes, Shrubs and TWL started by EBG in 7/2017; IRBR 10-11/18

# Converted to Python by IRBR 1/2019

# Version 5 includes:
#     + Dune diffusion technique
#     + Minimum dune elevation (mean high water line)



import numpy
import random
import math
import matplotlib.pyplot as plt
import time

Time = time.time()


############################################################################################################
# SET PARAMETERS       
############################################################################################################


################################
### TIME

TMAX = 10
StormStart = 3



################################
### COMPUTATIONAL DOMAIN (DECAMETERS)

DomainLength = 500
DomainWidth = 50
DuneDomain = numpy.zeros([TMAX, DomainLength])
ShrubDomainAll = numpy.zeros([TMAX, DomainLength, DomainWidth])
ShrubDomainFemale = numpy.zeros([TMAX, DomainLength, DomainWidth])
ShrubDomainMale = numpy.zeros([TMAX, DomainLength, DomainWidth])
ShrubPercentCover = numpy.zeros([TMAX, DomainLength, DomainWidth])



################################
### DUNE

# Initial height of the dune domain above (relative to) the NAVD88 MHW elevation on Hog
Dstart = 0.70
MHW = 0.46 # (m) Mean high water NAVD88 elevation; from EDM thesis; MHW + Dune Height = Dune Elevation
DuneDomain[0,:] = numpy.ones([1, DomainLength]) * (Dstart + (-0.1 + (0.1 - (-0.1)) * numpy.random.rand(1,DomainLength)))

# Maximum dune height
Dmaxel = 4 # Max dune elevation (m NAVD88)
Dmax = Dmaxel - MHW # Max dune height (m)

# Dune growth parameter
# Houser et al 0.05 ? r ? 0.55
rmin = 0.05
rmax = 0.55
growthparam = rmax + (rmax-rmin) * numpy.random.rand(1,DomainLength)



################################
### STORM

# Number Per Year
lambda_storm = 52.398 # For Weibull (temp)
k_storm = 4.5301   
numstorm = 50 # Optional, if you want all years to have same number

# Total Water Level
mu_storm = 1.1075 #1.6075 #1.1075
sigma_storm = 0.3353 #0.3353



################################
### SHRUB

# Shrub Growth from JCZ
# Minimum elevation of 2 m NAVD88 for shrub growth
Dshrub = 2 - MHW # Minimum dune height for shrub growth

# Percent cover change (year 1-10)
PC = [4, 8, 10, 15, 15, 20, 35, 80, 100, 100]

# Seeds produced per cell per year (Fecundity)
# Seedmin=1000; %???
# Seedmax=10000; % from Kwit et al 2004 Oecologia
Seedmin = 100 # Temp
Seedmax = 1000 # Temp

# Seed dispersal distance
# Dows 2016 has the data
# Took the 0 thicket data from 0-2m and then the grassland data from 50-300m. then fit a lognormal distribution (better than exponential) to data.
# Use fn 'lognrnd' to get dispersal distances. i.e., DispDist = lognrnd(mu,sigma)
mu = -0.721891
sigma = 1.5

# Germination rate - 60% from Young et al 1994
GermRate = 0.6

# Survival rate (from those reported in Lazarus and McGill)
# SurvRatemin = 10e-5
# SurvRatemax = 0.02
SurvRatemin = GermRate
SurvRatemax = GermRate

# Trees take 4-6 years to fruit; email from Julie
TimeFruit = 5

# Trees are 50/50 male/female - Hokkanen thesis 
Female = 0.5

# Initiate Shrubs
ShrubDomainFemale[0, round(DomainLength/2), round(DomainWidth/2)] = 5


############################################################################################################
# RUN MODEL      
############################################################################################################

for t in range(1, TMAX): # Yearly time steps - actual time = t + 1
    
    ### Print time step to screen
    ts = 'Time Step: '
    print("\r", 'Time Step: ', t+1, end = "")


    ### Dune Growth
    G = growthparam * DuneDomain[t-1,:] * (1 - DuneDomain[t-1,:] / Dmax)
    DuneDomain[t,:] = G + DuneDomain[t-1,:]
    
    ### Shrubs
    # Age the shrubs in years
    for k in range(DomainWidth): # Loop through each row of island width (i.e. from ocean to mainland side of island)
        # Find percent cover for each cell
        for j in range(DomainLength): # Loop through each shrub cell in domain
            ageF = ShrubDomainFemale[t-1,j,k] # Get age of shrub in cell from last time step
            ageM = ShrubDomainMale[t-1,j,k]
            if ageF > 0: # Shrub at this location is Female
                ShrubDomainFemale[t,j,k] = ShrubDomainFemale[t-1,j,k] + 1
                ShrubDomainMale[t,j,k] = 0
                if ageF > 10:
                   ShrubPercentCover[t,j,k] = 100
                else:
                   ShrubPercentCover[t,j,k] = PC[int(ageF)-1]
            elif ageM > 0: # Shrub at this location is Male
                ShrubDomainMale[t,j,k] = ShrubDomainMale[t-1,j,k] + 1
                ShrubDomainFemale[t,j,k] = 0
                if ageM > 10:
                    ShrubPercentCover[t,j,k] = 100
                else:
                    ShrubPercentCover[t,j,k] = PC[int(ageM)-1]
            else: # No shrub at this location
                ShrubDomainMale[t,j,k] = 0
                ShrubDomainFemale[t,j,k] = 0
                ShrubPercentCover[t,j,k] = 0
                
    
    ### Disperse seeds
    for k in range(DomainWidth): # Loop through each row of island width (i.e. from ocean to mainland side of island)
        # For all cells with a shrub
        FemaleShrubs = ShrubDomainFemale[:,:,k]
        I = [index for index, value in enumerate(FemaleShrubs[t, :]) if value >= TimeFruit]
        numShrubCells = len(I)
        # Determine how many seeds in each cell
        Seedspercell = numpy.random.randint(Seedmin,high = Seedmax, size = numShrubCells)
        # For each shrub cell, determine survival rate for the cell in this year
        SeedSurvYear = SurvRatemin + (SurvRatemax - SurvRatemin) * numpy.random.rand(numShrubCells) #<-- ,1)?

        for i in range(numShrubCells):
            # For each shrub cell producing seeds, generate a random # of random numbers each representing a single seed
            randseeds = numpy.random.rand(Seedspercell[i])
            # Find how many seeds produced in each cell actually survive
            Survivors = len(randseeds[randseeds < SeedSurvYear[i]]) #<--

            # Determine distance, rounding to nearest integer
            if Survivors > 0:
                DispDist = numpy.round(numpy.random.lognormal(mu,sigma,Survivors))
            else:
                DispDist = []

            # If there are seeds to disperse
            if len(DispDist) > 0:
                for j in range(len(DispDist+1)): # Loop through each individual seed to disperse
                    # Create a meshgrid to find coordinates of all points that are dropdistance from origin 
                    if DispDist[j] > 0:
                        gridsize = int(DispDist[j] * 2 + 2)
                        X, Y = numpy.meshgrid(range(1, gridsize), range(1, gridsize))
                        originX = I[i] # Sets coordinates of plant where seed is coming from
                        originY = k
                        matOriginX = math.ceil(len(X)/2) 
                        matOriginY = matOriginX
                        distMat = numpy.round(numpy.sqrt((X-matOriginX)**2 + (Y-matOriginY)**2)) # Find the distance from origin to every other point on island
                        # Find coordinates of each point on island that is dropdistance away from origin
                        coords = numpy.where(distMat == DispDist[j]) 
                        row = coords[0]
                        col = coords[1] 
                        
                        # Randomly selct one of those points as the target point - this means equal probability of dispersal in every direction
                        if len(col) == 0:
                            targetX = originX
                            targetY = originY
                        else:
                            xx = random.randint(0, (len(col))-1)
                            matTargetX = col[xx]
                            matTargetY = row[xx]
                            
                            targetX = originX + (matTargetX-matOriginX)
                            targetY = originY + (matTargetY-matOriginY)
                        
                        
                        # Put a plant in the ground if 
                        #   -the dropdistance>0, 
                        #   -the target drop location is within the island domain, 
                        #   -there is no shrub at the receiving cell (male or female),
                        #   -the receiving cell has a tall enough fronting dune
                        if  targetX > 0 and targetX < DomainLength and targetY > 0 and targetY < DomainWidth \
                            and ShrubDomainFemale[t,targetX,targetY] == 0 and ShrubDomainMale[t,targetX,targetY] == 0 \
                            and DuneDomain[t,targetX] >= Dshrub:
                            # Decide if the tree wll be a female or male
                            if random.random() > Female:
                                ShrubDomainFemale[t,targetX,targetY] = 1
                            else:
                                ShrubDomainMale[t,targetX,targetY] = 1
 
    
    
    ### Storms
    if t >= StormStart:
        # Select number of storms for this time step from weibull distribution
        numstorm = round(numpy.random.weibull(k_storm) * lambda_storm)
        
        if numstorm > 0:
            # Select TWL of each storm from normal distribution
            TWL = numpy.random.normal(mu_storm, sigma_storm, numstorm) # (m NAVD88) 
            
            # Overwash/Erosion
            # Reduce dune height for overwashed dunes & kill shrubs
            OWdist = numpy.random.randint(30, high=DomainWidth, size=DomainLength) # Random overwash penetration distance, a single distance for every storm of the year
            for n in range(numstorm):
                Low = numpy.random.uniform(0.3, 0.8) # Long et al. (2014) - 80% of observed change fall within 30-80% of pre-storm dune height
                Dow = [index for index, value in enumerate((DuneDomain[t, :] + MHW)) if value < TWL[n]]
                
                # Remove the relevant slice of the array to opertate on
                ShrubMortF = ShrubDomainFemale[t,:,:]
                ShrubMortM = ShrubDomainMale[t,:,:]
                for d in range(len(Dow)):
                    DuneDomain[t,Dow[d]] = DuneDomain[t,Dow[d]] * Low # Less efficient - but works
                    # If dune is lowered to essentially zero, allow for chance of regrowth by raising dune height to 5 cm
                    if DuneDomain[t,Dow[d]] < 0.05:
                        DuneDomain[t,Dow[d]] = 0.05
                    # Kill 1-year-old shrubs behind overwashed dunes
                    for w in range(OWdist[Dow[d]]):
                        if ShrubMortF[Dow[d], w] == 1 or ShrubMortM[Dow[d], w] == 1:
                            ShrubMortF[Dow[d],w] = 0
                            ShrubMortM[Dow[d],w] = 0
                # Put the relevant slice of the array back in
                ShrubDomainFemale[t,:,:] = ShrubMortF
                ShrubDomainMale[t,:,:] = ShrubMortM
                
                
                # Dune Height Diffusion - smoothes dune height by redistributing sand from high dune to neighboring low dune 
                AA = 0.8 # Maximum height difference between adjacent dune cells, adjustable
                for d in range(2,DomainLength): # Loop L to R
                    Ldiff = DuneDomain[t,d] - DuneDomain[t,d-1] # Calculate height difference
                    # Subtract sand from taller dune cell and add to adjacent shorter one (assumes sand within dunes is conserved)
                    if Ldiff > AA:
                        sub = (Ldiff - AA)/2 # Height difference in excess of maximum, divided by two
                        DuneDomain[t,d] = DuneDomain[t,d] - sub # Lower dune that's too tall
                        DuneDomain[t,d-1] = DuneDomain[t,d-1] + sub # Raise dune that's too short

                for d in range((DomainLength-2),0,-1): # Loop R to L
                    Rdiff = DuneDomain[t,d] - DuneDomain[t,d+1]
                    if Rdiff > AA:
                        sub = (Rdiff - AA)/2
                        DuneDomain[t,d] = DuneDomain[t,d] - sub
                        DuneDomain[t,d+1] = DuneDomain[t,d+1] + sub

        
        # Kill all exposed shrubs behind dunes lowered below threshold height
        for d in range(DomainLength):
            if DuneDomain[t,d] < Dshrub:
                ShrubDomainFemale[t,d,range(OWdist[d])] = 0
                ShrubDomainMale[t,d,range(OWdist[d])] = 0



ShrubDomainAll = ShrubDomainMale + ShrubDomainFemale

# End of model run
print()
print('Elapsed Time: ', time.time() - Time, 'sec') # Print elapsed time of simulation




############################################################################################################
# PLOT RESULTS      
############################################################################################################


# Dune Height Over Time
duneFig = plt.figure(figsize=(22,5))
ax = duneFig.add_subplot(111)
cax = ax.matshow(DuneDomain, origin='lower', cmap='bwr', aspect='auto')
ax.xaxis.set_ticks_position('bottom')
cbar = duneFig.colorbar(cax)
cbar.set_label('Dune Height Above MHW (m)', rotation=270)
plt.xlabel('Alongshore Distance (dm)')
plt.ylabel('Time (yr)')
plt.title('Dunes Over Time')
plt.show()

# Shrub Age Domain at Simulation End
ShrubAllRot90 = numpy.rot90(ShrubDomainAll[TMAX-1,:,:])
shrubFig1 = plt.figure(figsize=(22,3))
ax = shrubFig1.add_subplot(111)
cax = ax.matshow(ShrubAllRot90, origin='lower')
ax.xaxis.set_ticks_position('bottom')
cbar = shrubFig1.colorbar(cax)
cbar.set_label('Shrub Age (yr)', rotation=270)
plt.xlabel('Alongshore Distance (dm)')
plt.ylabel('Cross-Shore Diatance (dm)')
plt.title('Final Shrub Age')
plt.show()

## Percent Cover Domain at Simulation End
#ShrubPCRot90 = numpy.rot90(ShrubPercentCover[TMAX-1,:,:])
#shrubFig2 = plt.figure(figsize=(22,3))
#ax = shrubFig2.add_subplot(111)
#cax = ax.matshow(ShrubPCRot90, origin='lower', cmap='RdYlGn')
#ax.xaxis.set_ticks_position('bottom')
#cbar = shrubFig2.colorbar(cax)
#cbar.set_label('Shrub Percent Cover', rotation=270)
#plt.xlabel('Alongshore Distance (dm)')
#plt.ylabel('Cross-Shore Diatance (dm)')
#plt.title('Final Shrub Cover')
#plt.show()




























