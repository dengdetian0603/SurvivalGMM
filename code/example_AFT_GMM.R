source("./AFT_GMM.R")

beta = c(0.2, 0.5, -0.5)
sigma = 0.5
true.par = data.frame(value = c(beta, sigma),
                      par.name = c("b0", "b1", "b2", "sigma"))

sim.study.obj = SimStudyAFTGMM(n.sample = 100, tcut0 = 0.5, n.rep = 200,
                               max.core = 3, pr.cens = 0.75)

dt = as.data.table(sim.study.obj)
dt = dt[value < 5 & value > -5]
dt[, list(mean = mean(value), sd = sd(value)), by = method:par.name]

PlotSimStudy(dt, true.par)

save(dt, sim.study.obj,
     file = "../../Workspace/AFT_1tcut_n100_prcens0.75.RData")