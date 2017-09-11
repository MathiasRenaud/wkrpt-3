require(Kaphi)

# If working from GUI
setwd('~/git/MathiasRenaud/wkrpt-3')

# Coalescent Model
config <- load.config('YAML/coalescent.yaml')
config <- set.model(config, 'const.coalescent')

# Yule Model
#config <- load.config('YAML/yule.yaml')
#config <- set.model(config, 'yule')

# Birth-Death Model
#config <- load.config('YAML/bd.yaml')
#config <- set.model(config, 'bd')

set.seed(10)

obs.tree <- read.tree(file='cyano.nwk')
obs.tree <- parse.input.tree(obs.tree, config)

# initialize workspace
ws <- init.workspace(obs.tree, config)

# run ABC-SMC
# update file name and model!
res <- run.smc(ws, trace.file='cyano_coal_2.tsv', model='const.coalescent', verbose=TRUE)

trace <- read.table('cyano_yule_2.tsv', header=T, sep='\t')

par(mar=c(5,5,2,2))
plot(
  sapply(split(trace$lambda*trace$weight, trace$n), sum), 
  ylim=c(0, 2), 
  type='o',
  xlab='Iteration', 
  ylab='Mean lambda',
  cex.lab=1,
  main='Trajectory of Mean Lambda (Yule Model, 1000 particles)'
)
abline(h=0.16, lty=2)
