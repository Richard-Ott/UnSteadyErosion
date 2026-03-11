"""Interactive Python equivalent of Test_MCMC.m.

Edit the configuration block below and run this file directly.
The workflow mirrors the MATLAB script, but omits the plotting section.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat

from pyinversion import constants, forward_model, likelihood, mcmc_sampler, plots, sample_parameters, utils


def make_test_data(scenario: str, n: int) -> dict:
    """Port of make_test_data.m for the scenarios used in Test_MCMC.m."""
    testdata = {
        'lat': np.full(n, 30.0),
        'lon': np.full(n, 10.0),
        'altitude': np.full(n, 500.0),
    }

    if scenario == 'step':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([20, 50, 100, 300, 50, 400, 100,
                                  100, 2000, 4000, 200, 50, 100, 150], dtype=float)
        testdata['changeVariable'] = np.array([], dtype=float)
    elif scenario == 'samestep':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([20, 50, 100, 40, 60, 80, 90], dtype=float)
        testdata['chg'] = np.array([20.0])
        testdata['changeVariable'] = testdata['chg']
    elif scenario == 'samebackground_step':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([50.0])
        testdata['chg'] = np.array([10, 20, 5, 40, 2, 25, 30], dtype=float)
        testdata['changeVariable'] = testdata['chg']
    elif scenario == 'samebackground_samestep':
        testdata['t'] = np.array([1000.0])
        testdata['e'] = np.array([50.0])
        testdata['chg'] = np.array([10.0])
        testdata['changeVariable'] = testdata['chg']
    elif scenario == 'spike':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([20, 50, 100, 300, 30, 400, 50], dtype=float)
        testdata['loss'] = np.array([50, 20, 1, 30, 10, 20, 40], dtype=float)
        testdata['changeVariable'] = testdata['loss']
    elif scenario == 'samespike':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([20, 50, 100, 300, 50, 400, 50], dtype=float)
        testdata['loss'] = np.array([50.0])
        testdata['changeVariable'] = testdata['loss']
    elif scenario == 'samebackground_spike':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([50.0])
        testdata['loss'] = np.array([50, 10, 100, 30, 20, 50, 5], dtype=float)
        testdata['changeVariable'] = testdata['loss']
    elif scenario == 'samebackground_samespike':
        testdata['t'] = np.array([1500.0])
        testdata['e'] = np.array([50.0])
        testdata['loss'] = np.array([20.0])
        testdata['changeVariable'] = testdata['loss']
    elif scenario == 'curve':
        pollen = loadmat('data/pollen.mat', simplify_cells=True)['pollen']
        years_bp = np.asarray(pollen['yearsBP'], dtype=float)
        perc_tree = np.asarray(pollen['percTree'], dtype=float)

        timebreaks = np.array([10000, 6200, 700, 0], dtype=float)
        mean_perc_tree = []
        for i in range(len(timebreaks) - 1):
            mask = (years_bp >= timebreaks[i + 1]) & (years_bp < timebreaks[i])
            mean_perc_tree.append(np.mean(perc_tree[mask]))
        mean_perc_tree = np.asarray(mean_perc_tree, dtype=float)
        no_tree_perc = 100.0 - mean_perc_tree
        curvechanges = no_tree_perc[1:] / no_tree_perc[0] - 1.0

        testdata['t'] = timebreaks[1:-1]
        testdata['e'] = np.array([10, 500, 100, 50, 20, 80, 150], dtype=float)
        testdata['chg'] = np.array([10, 1, 4, 5, 0.1, 100, 20], dtype=float)
        testdata['curvechange'] = curvechanges
        testdata['changeVariable'] = testdata['chg']
    else:
        raise ValueError(f'Unknown scenario: {scenario}')

    testdata['steps'] = len(testdata['t'])
    return testdata


def build_init_params(prior_range: np.ndarray, var_names: list[str], num_chains: int, seed: int) -> dict:
    """Sample initial chain positions from the prior ranges."""
    rng = np.random.default_rng(seed)
    init_params = {}
    for index, name in enumerate(var_names):
        init_params[name] = rng.uniform(prior_range[index, 0], prior_range[index, 1], size=(num_chains,))
    return init_params


def reshape_samples(samples: dict, var_names: list[str], nwalks: int, num_samples: int) -> np.ndarray:
    """Convert merged NumPyro samples dict into MATLAB-like stacked array."""
    return np.stack(
        [np.asarray(samples[name]).reshape((nwalks, num_samples)) for name in var_names],
        axis=0,
    )


# -----------------------------------------------------------------------------
# Configuration block: edit these values and run the script interactively.
# -----------------------------------------------------------------------------
scenario = 'step'
n = 7
nWalks = 3
nsamples = 10000
burnin = 0.2
seed = 0
export = False
filetag = 'test'
make_plots = True
show_plots = True

tdata = make_test_data(scenario, n)
nlogical = np.column_stack([
    np.ones(n, dtype=bool),
    np.ones(n, dtype=bool),
    np.zeros(n, dtype=bool),
])

# Priors 
t_prior = [1e2, 6e3]
e_prior = [0, 4.5e3]
loss_prior = [0, 200]
chg_prior = [0, 50]
prior_range, var_names = likelihood.make_prior_and_varnames(
    scenario, t_prior, e_prior, loss_prior, chg_prior, n, tdata['steps']
)

consts_obj = constants.load_constants('online-calculators-v3/consts_v3.mat')
sp = sample_parameters.sample_parameters(tdata['lat'], tdata['lon'], tdata['altitude'], consts_obj)
if scenario == 'curve':
    sp['curvechange'] = np.asarray(tdata['curvechange'], dtype=float)
    sp['t'] = np.asarray(tdata['t'], dtype=float)

sp_jax = {key: jnp.asarray(value) if isinstance(value, np.ndarray) else value for key, value in sp.items()}
consts_jax = {key: value for key, value in consts_obj.items()}


def measured_nuclide_indices(nlogical_array: np.ndarray) -> np.ndarray:
    n_samp = nlogical_array.shape[0]
    return np.concatenate([
        np.where(nlogical_array[:, 0])[0],
        n_samp + np.where(nlogical_array[:, 1])[0],
        2 * n_samp + np.where(nlogical_array[:, 2])[0],
    ])


use_jax_step_forward = scenario == 'step'
measured_indices = measured_nuclide_indices(nlogical)


def forward_model_fn(model_vector):
    if use_jax_step_forward:
        full_prediction = forward_model.Nforward_wrapper_jax(model_vector, sp_jax, consts_jax, scenario, tdata['steps'], nlogical)
        return full_prediction[measured_indices]
    return forward_model.Nforward_wrapper(model_vector, sp_jax, consts_jax, scenario, tdata['steps'], nlogical)


forward_model_fn = jax.jit(forward_model_fn)


mtest = np.concatenate([
    np.asarray(tdata['t'], dtype=float).ravel(),
    np.asarray(tdata['e'], dtype=float).ravel(),
    np.asarray(tdata['changeVariable'], dtype=float).ravel(),
])
test_obs = np.asarray(forward_model_fn(jnp.asarray(mtest)))
dtest_obs = np.maximum(np.abs(test_obs) * 0.08, 1e-12)

model_fn = likelihood.create_numpyro_model(
    jnp.asarray(test_obs),
    jnp.asarray(dtest_obs),
    forward_model_fn,
    nlogical,
    prior_range,
    var_names,
    use_safe_jacfwd=not use_jax_step_forward,
)

num_warmup = int(burnin * nsamples)
num_samples = max(1, nsamples - num_warmup)
init_params = build_init_params(prior_range, var_names, nWalks, seed)

rng_key = jax.random.PRNGKey(seed)
samples, sampler = mcmc_sampler.run_mcmc_inversion(
    model_fn,
    num_chains=nWalks,
    num_samples=num_samples,
    num_warmup=num_warmup,
    rng_key=rng_key,
    init_params=init_params,
)

samples_array = reshape_samples(samples, var_names, nWalks, num_samples)

flat_models = jnp.asarray(samples_array.transpose(1, 2, 0).reshape((-1, samples_array.shape[0])))
batched_forward_model_fn = jax.jit(jax.vmap(forward_model_fn))
flat_predictions = np.asarray(batched_forward_model_fn(flat_models))
flat_log_like = np.sum(
    -0.5 * ((test_obs[None, :] - flat_predictions) / dtest_obs[None, :]) ** 2
    - np.log(np.sqrt(2 * np.pi) * dtest_obs[None, :]),
    axis=1,
)
log_like = flat_log_like.reshape((nWalks, num_samples))

best_chain, best_sample = np.unravel_index(np.argmax(log_like), log_like.shape)
best_model = samples_array[:, best_chain, best_sample]
best_pred = np.asarray(forward_model_fn(jnp.asarray(best_model)))

figures = {}
if make_plots:
    figures['autocorrelation'], autocorr = plots.autocorrelationplot(samples_array)
    figures['chains'] = plots.chainplot(samples_array, var_names, prior_range, true_vals=mtest)
    figures['cornerplot'] = plots.ecornerplot(
        samples_array,
        ks=True,
        color=(0.3, 0.3, 0.3),
        names=var_names,
        bestmodel=best_model,
        truevals=mtest,
    )
    figures['barplot'] = plots.barplot_parameters(
        samples_array,
        var_names,
        prior_range,
        bestmodel=best_model,
        truevals=mtest,
    )
    n10_obs, n14_obs, _ = plots.split_by_nuclide(test_obs, nlogical)
    d10_obs, d14_obs, _ = plots.split_by_nuclide(dtest_obs, nlogical)
    n10_pred, n14_pred, _ = plots.split_by_nuclide(best_pred, nlogical)
    best_pred_plot = np.concatenate([n10_pred, n14_pred])
    figures['datafit'] = plots.conc_modelledVSobserved(best_pred_plot, n10_obs, d10_obs, n14_obs, d14_obs)

print('-' * 70)
print(f'Best chain: {best_chain}')
print(f'Best sample: {best_sample}')
print(f'Best log-likelihood: {log_like[best_chain, best_sample]:.4f}')
print(f'Mean absolute residual: {np.mean(np.abs(best_pred - test_obs)):.6g}')

if export:
    output_dir = Path('output')
    output_dir.mkdir(exist_ok=True)
    stem = f'{filetag}_{scenario}_python'
    mat_path = output_dir / f'{stem}.mat'
    pkl_path = output_dir / f'{stem}.pkl'

    savemat(
        mat_path,
        {
            'scenario': np.array([scenario], dtype=object),
            'mtest': mtest,
            'testObs': test_obs,
            'testObs_sigma': dtest_obs,
            'prior_range': prior_range,
            'var_names': np.array(var_names, dtype=object),
            'samples_array': samples_array,
            'logLike': log_like,
            'best_model': best_model,
            'best_pred': best_pred,
            'Nlogical': nlogical,
        },
    )

    with open(pkl_path, 'wb') as handle:
        pickle.dump(
            {
                'scenario': scenario,
                'tdata': tdata,
                'mtest': mtest,
                'testObs': test_obs,
                'testObs_sigma': dtest_obs,
                'prior_range': prior_range,
                'var_names': var_names,
                'samples_array': samples_array,
                'logLike': log_like,
                'best_model': best_model,
                'best_pred': best_pred,
                'Nlogical': nlogical,
            },
            handle,
        )
    if figures:
        plots.save_diagnostic_figures(figures, output_dir, f'{filetag}_{scenario}')
    print(f'Saved results to {mat_path} and {pkl_path}')

if figures and show_plots:
    plt.show()