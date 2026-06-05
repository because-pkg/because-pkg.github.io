import jax
import jax.numpy as jnp

def test_exec():
    def model(x):
        local_env = {"jnp": jnp, "x": x}
        code = """
y = jnp.sin(x) * 2.0
"""
        exec(code, globals(), local_env)
        return local_env["y"]
    
    jitted = jax.jit(model)
    print("Result:", jitted(1.0))

test_exec()
