# Inference benchmark (`infbench`)

The goal of `infbench` is to compare various sample-efficient approximate inference algorithms which have been proposed in the machine learning literature so as to deal with (moderately) expensive likelihoods. In particular, we want to infer the posterior over model parameters and (an approximation of) the model evidence or marginal likelihood, that is the normalization constant of the posterior. Crucially, we assume that the budget of likelihood evaluations is of the order of a few hundreds.

Notably, this goal is more ambitious than simply finding the maximum of the posterior (MAP), problem that we previously tackled with [Bayesian Adaptive Direct Search](https://github.com/lacerbi/bads) (aka BADS).

Our preliminary benchmark - to be released soon - shows that existing inference algorithms perform quite poorly at reconstructing posteriors (or evaluating their normalization constant) with both syntethic and real pdfs that have moderately challenging features, showing that this is a much harder problem.

This project is currently work in progress, stay tuned for updates or write to me at <luigi.acerbi@unige.ch>.

## How to run the benchmark

You can run the benchmark on one test problem as follows:
```
> options = struct('SpeedTests',false);
> [probstruct,history] = infbench_run('vbmc18',testfun,D,[],algo,id,options);
```
The arguments are:

- `testfun` (string) is the test pdf, which can be `'cigar'`, `'lumpy'`, `'studentt'` (synthetic test functions with different properties), or `'goris2015'` (real model fitting problem with neuronal data, see [here](https://github.com/lacerbi/infbench/tree/master/matlab/problems/goris2015)).
- `D` (integer) is the dimensionality of the problem for the synthetic test functions (typical values are from `D=2` to `D=10`), and `7` or `8` for the `goris2015` test set (corresponding to two different neuron datasets).
- `algo` (string) is the inference algorithm being tested. Specific settings for the chosen inference algorithm are selected with `'algo@setting'`, e.g. `'agp@debug'`.
- `id` (integer) is the run id, used to set the random seed. It can be a vector for multiple consecutive runs.
- `options` (struct) sets various options for the benchmark. For a fast test, I recommend to set the field `SpeedTests` to `false` since initial speed tests can be quite time consuming.

And the outputs are:

- `probstruct` (struct) describes the current problem in detail.
- `history` (struct array) summary statistics of the run for each id.
