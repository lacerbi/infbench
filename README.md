# Introduction

This repository contains code for the inference benchmarks and results reported in several papers on *Variational Bayesian Monte Carlo* (VBMC) [[1,2,3](#references)]. 

The Variational Bayesian Monte Carlo (VBMC) algorithm is available at [this repository](https://github.com/lacerbi/vbmc).

# Inference benchmark (`infbench`)

The goal of `infbench` is to compare various sample-efficient approximate inference algorithms which have been proposed in the machine learning literature so as to deal with (moderately) expensive and potentially noisy likelihoods. In particular, we want to infer the posterior over model parameters and (an approximation of) the model evidence or marginal likelihood, that is the normalization constant of the posterior. Crucially, we assume that the budget is of the order of up to several hundreds likelihood evaluations.

Notably, this goal is more ambitious than simply finding the maximum of the posterior (MAP), problem that we previously tackled with [Bayesian Adaptive Direct Search](https://github.com/lacerbi/bads) (aka BADS).

Our [first benchmark](#how-to-run-the-original-benchmark-vbmc18) shows that existing inference algorithms perform quite poorly at reconstructing posteriors (or evaluating their normalization constant) with both syntethic and real pdfs that have moderately challenging features, showing that this is a much harder problem [[1,2](#reference)].

Our [second extensive benchmark](#how-to-run-the-extensive-benchmark-vbmc20) shows that the latest version of VBMC (v1.0, June 2020) beats other state-of-the-art methods on real computational problems when dealing with *noisy* log-likelihood evaluations, such as those arising from simulation-based estimation techniques [[3](#reference)].

## How to run the original benchmark (`vbmc18`)

You can run the benchmark on one test problem in the `vbmc18` problem set as follows:
```
> options = struct('SpeedTests',false);
> [probstruct,history] = infbench_run('vbmc18',testfun,D,[],algo,id,options);
```
The arguments are:

- `testfun` (string) is the test pdf, which can be `'cigar'`, `'lumpy'`, `'studentt'` (synthetic test functions with different properties), or `'goris2015'` (real model fitting problem with neuronal data, see [here](https://github.com/lacerbi/infbench/tree/master/matlab/problems/goris2015)).
- `D` (integer) is the dimensionality of the problem for the synthetic test functions (typical values are from `D=2` to `D=10`), and `7` or `8` for the `goris2015` test set (corresponding to two different neuron datasets).
- `algo` (string) is the inference algorithm being tested. Specific settings for the chosen inference algorithm are selected with `'algo@setting'`, e.g. `'agp@debug'`.
- `id` (integer) is the run id, used to set the random seed. It can be an array, such as `id=1:5`, for multiple consecutive runs.
- `options` (struct) sets various options for the benchmark. For a fast test, I recommend to set the field `SpeedTests` to `false` since initial speed tests can be quite time consuming.

The outputs are:

- `probstruct` (struct) describes the current problem in detail.
- `history` (struct array) summary statistics of the run for each id.

## How to run the extensive benchmark (`vbmc20`)

The `vbm20` benchmark includes a number of real, challenging models and data largely from computational and cognitive neuroscience, from `D = 3` to `D = 9`. The benchmark is mostly designed to test methods that deal with *noisy* log-likelihood evaluations.

You can run the benchmark on one test problem in the `vbmc20` problem set as follows:
```
> options = struct('SpeedTests',false);
> [probstruct,history] = infbench_run('vbmc20',testfun,D,[],algo,id,options);
```
The arguments are:

- `testfun` (string) indicates the tested model, which can be `'wood2010'` (Ricker), `'krajbich2010'` (aDDM), `'acerbi2012'` (Timing), `'acerbidokka2018'` (Multisensory), `'goris2015b'` (Neuronal), `'akrami2018b'` (Rodent). The additional model presented in the Supplement is `'price2018'` (g-and-k).
- `D` (integer) is the dataset, with `D = 1` for Ricker, Timing, Rodent, g-and-k; `D = 1` or `D = 2` for aDDM and Multisensory; `D = 7` or `D = 8` for Neuronal.
- `algo` (string) is the inference algorithm being tested. For the algorithms tested in the paper [[3](#references)], use `'vbmc@imiqr'`, `'vbmc@viqr'`, `'vbmc@npro'`, `'vbmc@eig'`, `'parallelgp@v3'`, `'wsabiplus@ldet'`.

For the other input and output arguments, see [above](#how-to-run-the-original-benchmark-vbmc18).

Code used to generate figures in the paper [[3](#references)] is available in [this folder](https://github.com/lacerbi/infbench/tree/master/matlab/figs/vbmc_paper2020). However, you will first need to run the benchmark (due to space limitations, we cannot upload the bulk of the numerical results here).

## References

1. Acerbi, L. (2018). Variational Bayesian Monte Carlo. In *Advances in Neural Information Processing Systems 31*: 8222-8232. ([paper + supplement on arXiv](https://arxiv.org/abs/1810.05558), [NeurIPS Proceedings](https://papers.nips.cc/paper/8043-variational-bayesian-monte-carlo))
2. Acerbi, L. (2019). An Exploration of Acquisition and Mean Functions in Variational Bayesian Monte Carlo. In *Proc. Machine Learning Research* 96: 1-10. 1st Symposium on Advances in Approximate Bayesian Inference, Montréal, Canada. ([paper in PMLR](http://proceedings.mlr.press/v96/acerbi19a.html))
3. Acerbi, L. (2020). Variational Bayesian Monte Carlo with Noisy Likelihoods. arXiv preprint (to appear).

### Contact

This repository is currently actively used for research, stay tuned for updates:

- [Follow me on Twitter](https://twitter.com/AcerbiLuigi) for updates about my work on model inference and other projects I am involved with;
- If you have questions or comments about this work, get in touch at <luigi.acerbi@unige.ch> (putting 'infbench' in the subject of the email).

### License

The inference benchmark is released under the terms of the [GNU General Public License v3.0](https://github.com/lacerbi/infbench/blob/master/LICENSE.txt).
