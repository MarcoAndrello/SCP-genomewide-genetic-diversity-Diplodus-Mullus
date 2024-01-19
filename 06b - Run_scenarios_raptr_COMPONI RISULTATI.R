prob_50gs_up <- prob_50gs
res_50gs_up <- res_50gs
prob_50gs_up[1:5] <- problems_raptr[1:5]
res_50gs_up[1:5] <- results_raptr[1:5]
prob_50gs_up -> prob_50gs
res_50gs_up -> res_50gs
prob_50gs
res_50gs
save(prob_50gs, res_50gs, file="Results/Results_raptr_50gs.RData")

rm(list=ls())
prob_20gs_up <- prob_20gs
res_20gs_up <- res_20gs
prob_20gs_up[6:12] <- prob_20gs[6:12]
res_20gs_up[6:12] <- res_20gs[6:12]
prob_20gs_up[1:5] <- problems_raptr[6:10]
res_20gs_up[1:5] <- results_raptr[6:10]
prob_20gs_up -> prob_20gs
res_20gs_up -> res_20gs
prob_20gs
res_20gs
save(prob_20gs, res_20gs, file="Results/Results_raptr_20gs.RData")

### SE TUTTO OK, CANCELLARE:
# Results_raptr.RData, Results_raptr_50gs_perm6.RData etc
## MI SERVONO I PROBLEMS IN QUESTO FILE DI RISULTATI??