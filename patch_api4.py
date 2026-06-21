import sys
with open("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/api.py", "r") as f:
    content = f.read()

replacement = """        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*There are not enough devices.*")
            test_mcmc = MCMC(kernel_base, num_warmup=num_warmup, num_samples=num_samples, num_chains=dsep_chains, thinning=thinning, progress_bar=False)
            test_mcmc.run(subkey, **jax_dsep_data)"""

content = content.replace("""        test_mcmc = MCMC(kernel_base, num_warmup=num_warmup, num_samples=num_samples, num_chains=dsep_chains, thinning=thinning, progress_bar=False)
        
        test_mcmc.run(subkey, **jax_dsep_data)""", replacement)

with open("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/api.py", "w") as f:
    f.write(content)
