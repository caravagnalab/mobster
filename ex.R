# Example [mobster_0.1.0.tar.gz]

# It requires you to install first
# remotes::install_github("caravagn/pio")
# remotes::install_github("caravagn/easypar")

require(mobster)


pi = c(`Tail` = .3, `C1` = .2, `C2` = .5)
N = round(pi * 500)

tail_samples = sads::rpareto(N['Tail'], shape = 2, scale = 0.05)
C1_samples = rbinom(N['C1'], size = 100, prob = .3)/100
C2_samples = rbinom(N['C2'], size = 100, prob = .5)/100

input = data.frame(
  VAF = c(tail_samples, C1_samples, C2_samples),
  label = c(
    rep('Tail', N['Tail']),
    rep('C1', N['C1']),
    rep('C2', N['C2'])
  )
)

# Input data
ggplot(input, aes(VAF, fill = label)) + geom_histogram(binwidth = 0.01) + xlim(0,1)

# Fit
fit = mobster_fit(input, K = 1:3, parallel = TRUE, samples = 1, tail = c(TRUE, F))

# S3 methods
fit$best
plot(fit$best)

# Other plots
plot_latent_variables(fit$best)
plot_mixing_proportions(fit$best)

# other possible fits
fig = plot_model_selection(fit)
ggsave('example.pdf', fig, width = 10, height = 12)

# Clustering assignment(s)
fit$best$data

plot.dbpmm(fit$best, cutoff_assignment = .9)
ggplot(mobster:::Clusters(fit$best, cutoff_assignment = .99), aes(VAF, fill = label)) + geom_histogram(binwidth = 0.01)

# 10 resamples
bootstrap_results = mobster_bootstrap(fit_example$best, n.resamples = 600)

# Statistics ...
bootstrap_statistics = bootstrapped_statistics(fit_example$best, bootstrap_results = bootstrap_results)

boot_results$bootstrap_model
boot_results$bootstrap_statistics
boot_results$bootstrap_co_clustering

plot_bootstrap_model_frequency(bootstrap_results = bootstrap_results, bootstrap_statistics = boot_results)  
plot_bootstrap_mixing_proportions(bootstrap_results = bootstrap_results, bootstrap_statistics = boot_results)  


pdf("Bootstraps.pdf", width = 3, height = 2.5)
plot(fit_example$best)
lapply(boot$fits, plot)
dev.off()


# 10 resamples
x = random_dataset(seed = 123, Beta_variance_scaling = 100, N = 500)
x = mobster_fit(x$data, epsilon = 1e-5)
plot(x$best)
bootstrap_results = mobster_bootstrap(x$best, n.resamples = 15, epsilon = 1e-5)

bootstrap_results$fits = bootstrap_results$fits[!sapply(bootstrap_results$fits, function(e) inherits(e, 'error'))]

# Statistics ...
bootstrap_statistics = bootstrapped_statistics(x$best, bootstrap_results = bootstrap_results)

bootstrap_statistics$bootstrap_model
bootstrap_statistics$bootstrap_statistics
bootstrap_statistics$bootstrap_co_clustering

plot_bootstrap_model_frequency(bootstrap_results = bootstrap_results, bootstrap_statistics = bootstrap_statistics)  
plot_bootstrap_mixing_proportions(x$best, bootstrap_results = bootstrap_results, bootstrap_statistics = bootstrap_statistics)  
plot_bootstrap_Beta(x$best, bootstrap_results = bootstrap_results, bootstrap_statistics = bootstrap_statistics)  
plot_bootstrap_coclustering(x$best, bootstrap_results = bootstrap_results, bootstrap_statistics = bootstrap_statistics)  


pdf("Bootstraps.pdf", width = 3, height = 2.5)
plot(x$best)
lapply(bootstrap_results$fits, plot)
dev.off()




