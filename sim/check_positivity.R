# 2b) Checking for potential positivity issues by age
plot(density(filter(cens_data, cens==0 & atwork==1)$age), col="black", main="Age of active work", xlab="Age", ylab="Kernel Density", lwd=2)
lines(density(filter(sim_cohort, atwork==1)$age), col="red", lwd=2)
legend("topright", legend=c("Uncensored (limit=2)", "Cohort"), col=c("black", "red"), lty=1, lwd=2)

# also check by year
plot(density(filter(cens_data, cens==0 & atwork==1)$year), col="black", main="Year of active work", xlab="Year", ylab="Kernel Density", lwd=2)
lines(density(filter(sim_cohort, atwork==1)$year), col="red", lwd=2)
legend("topright", legend=c("Uncensored (limit=2)", "Cohort"), col=c("black", "red"), lty=1, lwd=2)