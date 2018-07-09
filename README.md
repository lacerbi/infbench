# Inference benchmark (`infbench`)

Benchmark of posterior and model inference algorithms for (moderately) expensive likelihoods.

The goal of `infbench` is to compare various sample-efficient approximate inference algorithms which have been proposed in the machine learning literature so as to deal with (moderately) expensive likelihoods. In particular, we want to infer the posterior over model parameters and (an approximation of) the model evidence or marginal likelihood. The basic assumption is that the number of available likelihood evaluations is of the order of a few hundreds.

Unfortunately, most of the tested algorithms perform quite poorly on synthetic posterior pdfs with challenging but realistic features, showing that this is a hard problem.

This project is currently work in progress, stay tuned for updates or write to me at <luigi.acerbi@unige.ch>.
