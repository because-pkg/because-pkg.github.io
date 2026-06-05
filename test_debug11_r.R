library(reticulate)
py_run_string("
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

def phylo_transform(numpyro, jnp, jax, dist, group_name, num_groups, matrix, z_raw, sigma):
    lambda_val = numpyro.sample(f'lambda_{group_name}', dist.Uniform(0, 1))
    I = jnp.eye(num_groups)
    V = lambda_val * matrix + (1.0 - lambda_val) * I
    L = jax.scipy.linalg.cholesky(V, lower=True)
    z_group = jnp.dot(L, z_raw)
    numpyro.deterministic(f'sigma_phylo_{group_name}', sigma * jnp.sqrt(lambda_val))
    return z_group

def model():
    z_raw = numpyro.sample('z_raw', dist.Normal(0, 1).expand([3]))
    sigma = numpyro.sample('sigma', dist.HalfNormal(1))
    matrix = jnp.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.5], [0.2, 0.5, 1.0]])
    z = phylo_transform(numpyro, jnp, jax, dist, 'phylo_DD', 3, matrix, z_raw, sigma)
    numpyro.sample('obs', dist.Normal(z, 0.1), obs=jnp.array([1.0, 2.0, 3.0]))

mcmc = MCMC(NUTS(model), num_warmup=10, num_samples=10)
mcmc.run(jax.random.PRNGKey(0))
samples = mcmc.get_samples()
print(samples.keys())
")
