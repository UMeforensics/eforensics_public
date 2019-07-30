### Michigan County Data ###

michigan_county <- read.csv("test/countypres_2000-2016.csv");  
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
summary(samples); 
