library(eforensics);

### PART 1: Running the Current Vignette 
## model
## -----
model    = 'bl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=7000, nCov=1, model=model)
data     = sim_data$data

## mcmc parameters
## ---------------
mcmc    = list(burnin=1, n.adapt=100, n.iter=100, n.chains=2)

## samples
## -------
help(eforensics)
message(samples    = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc))

message(summary(samples))
message(summary(samples, join.chains=T))

### Notice that this throws an error: ``Error in check_mcmc(mcmc)'' 
### This is from the fact that ``burnin'' should probably be ``burn.in''. 


### PART 2: Vignette, Corrected
## model
## -----
model    = 'bl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=7000, nCov=1, model=model)
data     = sim_data$data

## mcmc parameters
## ---------------
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2)

## samples
## -------
help(eforensics)
samples    = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc);

summary(samples); 
summary(samples, join.chains=T); 


### PART 3: Help Files
help(ef_simulateData);
help(eforensics); 


### PART 4: Generating Own Data
model = 'bl'
sim_data = ef_simulateData(n=2000, nCov=3, model=model)
data = sim_data$data
mcmc = list(burn.in=100, n.adapt=100, n.iter=1000, n.chains=2)
samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc)


### PART 5: Missing Covariates
model = 'bl'
sim_data = ef_simulateData(n=2000, nCov=1, model=model);
data = sim_data$data;
missing_rows <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$x1[missing_rows] <- NA;
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
message(samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc));


### PART 6: Missing Winning Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$w[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
message(samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc));


### PART 7: Missing Abstention Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$a[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
message(samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc));


### PART 8: Missing Total Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$N[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc);


### PART 9: Michigan County Data
michigan_county <- read.csv("countypres_2000-2016.csv");  
michigan_county <- michigan_county[which(michigan_county$year==2016 & michigan_county$state_po=="MI"),]
drop_vars <- names(michigan_county) %in% c("office", "FIPS", "year", "state", "state_po", "version"); 
michigan_county <- michigan_county[!drop_vars]; 
michigan_county_reshaped <- reshape(michigan_county, varying = NULL, idvar = "county", timevar = "candidate", direction = "wide"); 
michigan_county_data_reshaped_vars <- names(michigan_county_reshaped) %in% c("county", "candidatevotes.Donald Trump", "totalvotes.Other"); 
michigan_county_final <- michigan_county_reshaped[michigan_county_data_reshaped_vars];
names(michigan_county_final) <- c("county", "trump", "total"); 
county_registered_voters <- read.table("VoterCount.txt", header = FALSE, dec = "."); 
michigan_county_final$county <- gsub('\\s+', '', michigan_county_final$county)
names(county_registered_voters) <- c("county", "total_registered");
michigan_county_all <- merge(michigan_county_final, county_registered_voters, by = "county"); michigan_county_all$a <- michigan_county_all$total_registered - michigan_county_all$total;

names(michigan_county_all) <- c("county", "w", "total_votes", "N", "a"); 

mcmc = list(burn.in=10, n.adapt=100, n.iter=1000, n.chains=2);
samples = eforensics(w ~ N, a ~ N, data=michigan_county_all, model='bl', mcmc=mcmc);


### PART 10: King County, Washington Data (FAILS)
kingcounty_raw <- read.csv("kingcounty_data.csv"); 
kingcounty_pres <- kingcounty_raw[which(kingcounty_raw$Race=="US President & Vice President"),]; 
vars_to_drop <- names(kingcounty_pres) %in% c("LEG", "CC", "CG", "Race", "CounterGroup", "Party"); 
kingcounty_pres <- kingcounty_pres[!vars_to_drop]; 
kingcounty_pres_reshaped <- reshape(kingcounty_pres, varying = NULL, idvar = "Precinct", timevar = "CounterType", direction = "wide"); 
reshaped_vars <- names(kingcounty_pres_reshaped) %in% c("Precinct", "SumOfCount.Registered Voters", "SumOfCount.Donald J. Trump & Michael R. Pence"); 
kingcounty_pres_final <- kingcounty_pres_reshaped[reshaped_vars];
names(kingcounty_pres_final) <- c("precinct", "N", "trump"); 
kingcounty_pres_final$a <- kingcounty_pres_final$N - kingcounty_pres_reshaped$`SumOfCount.Darrell L. Castle & Scott N. Bradley` - kingcounty_pres_reshaped$`SumOfCount.Hillary Clinton & Tim Kaine` - kingcounty_pres_reshaped$`SumOfCount.Jill Stein & Ajamu Baraka` - kingcounty_pres_reshaped$`SumOfCount.Gary Johnson & Bill Weld` - kingcounty_pres_reshaped$`SumOfCount.Gloria Estela La Riva & Eugene Puryear` - kingcounty_pres_reshaped$`SumOfCount.Alyson Kennedy & Osborne Hart` - kingcounty_pres_reshaped$`SumOfCount.Times Blank Voted` - kingcounty_pres_reshaped$`SumOfCount.Donald J. Trump & Michael R. Pence` - kingcounty_pres_reshaped$`SumOfCount.Write-In`; 

mcmc = list(burn.in=10, n.adapt=100, n.iter=1000, n.chains=2);
message(samples = eforensics(trump ~ N, a ~ N, data=kingcounty_pres_final, model='bl', mcmc=mcmc));


### PART 11: King County, Washington Data Dropping 0's 
updated_kingcounty_pres_final <- kingcounty_pres_final[which(kingcounty_pres_final$N != 0),]; 
mcmc = list(burn.in=10, n.adapt=100, n.iter=1000, n.chains=2);
samples = eforensics(trump ~ N, a ~ N, data=updated_kingcounty_pres_final, model='bl', mcmc=mcmc);
