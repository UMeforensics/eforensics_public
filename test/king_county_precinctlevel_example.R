### sKing County, Washington Data (FAILS)

kingcounty_raw <- read.csv("test/kingcounty_data.csv"); 
kingcounty_pres <- kingcounty_raw[which(kingcounty_raw$Race=="US President & Vice President"),]; 
vars_to_drop <- names(kingcounty_pres) %in% c("LEG", "CC", "CG", "Race", "CounterGroup", "Party"); 
kingcounty_pres <- kingcounty_pres[!vars_to_drop]; 
kingcounty_pres_reshaped <- reshape(kingcounty_pres, varying = NULL, idvar = "Precinct", timevar = "CounterType", direction = "wide"); 
reshaped_vars <- names(kingcounty_pres_reshaped) %in% c("Precinct", "SumOfCount.Registered Voters", "SumOfCount.Donald J. Trump & Michael R. Pence"); 
kingcounty_pres_final <- kingcounty_pres_reshaped[reshaped_vars];
names(kingcounty_pres_final) <- c("precinct", "N", "trump"); 
kingcounty_pres_final$a <- kingcounty_pres_final$N - kingcounty_pres_reshaped$`SumOfCount.Darrell L. Castle & Scott N. Bradley` - kingcounty_pres_reshaped$`SumOfCount.Hillary Clinton & Tim Kaine` - kingcounty_pres_reshaped$`SumOfCount.Jill Stein & Ajamu Baraka` - kingcounty_pres_reshaped$`SumOfCount.Gary Johnson & Bill Weld` - kingcounty_pres_reshaped$`SumOfCount.Gloria Estela La Riva & Eugene Puryear` - kingcounty_pres_reshaped$`SumOfCount.Alyson Kennedy & Osborne Hart` - kingcounty_pres_reshaped$`SumOfCount.Times Blank Voted` - kingcounty_pres_reshaped$`SumOfCount.Donald J. Trump & Michael R. Pence` - kingcounty_pres_reshaped$`SumOfCount.Write-In`; 

mcmc = list(burn.in=10, n.adapt=100, n.iter=1000, n.chains=2);
print(try(eforensics(trump ~ N, a ~ N, data=kingcounty_pres_final, model='bl', mcmc=mcmc)));


### King County, Washington Data Dropping 0's 
updated_kingcounty_pres_final <- kingcounty_pres_final[which(kingcounty_pres_final$N != 0),]; 
mcmc = list(burn.in=10, n.adapt=100, n.iter=1000, n.chains=2);
samples = eforensics(trump ~ N, a ~ N, data=updated_kingcounty_pres_final, model='bl', mcmc=mcmc);
print(summary(samples));