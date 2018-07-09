# Inference benchmark (`infbench`)

The goal of `infbench` is to compare various sample-efficient approximate inference algorithms which have been proposed in the machine learning literature so as to deal with (moderately) expensive likelihoods. In particular, we want to infer the posterior over model parameters and (an approximation of) the model evidence or marginal likelihood, that is the normalization constant of the posterior. Crucially, we assume that the budget of likelihood evaluations is of the order of a few hundreds.

Notably, this goal is more ambitious than simply finding the maximum of the posterior (MAP), problem that we previously tackled with [Bayesian Adaptive Direct Search](https://github.com/lacerbi/bads) (aka BADS).

Our preliminary benchmark - to be released soon - shows that existing inference algorithms perform quite poorly at reconstructing posteriors (or evaluating their normalization constant) with both syntethic and real pdfs that have moderately challenging features, showing that this is a much harder problem.

This project is currently work in progress, stay tuned for updates or write to me at <luigi.acerbi@unige.ch>.
