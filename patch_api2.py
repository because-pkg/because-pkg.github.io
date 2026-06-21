import sys
with open("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/api.py", "r") as f:
    content = f.read()

new_content = content.replace(
    'print(f"dsep_equations_to_run: {dsep_equations_to_run}")\n    if dsep_equations_to_run is not None and len(dsep_equations_to_run) > 0 and isinstance(dsep_equations_to_run[0], dict):',
    'print(f"Type: {type(dsep_equations_to_run[0])}")\n    if dsep_equations_to_run is not None and len(dsep_equations_to_run) > 0 and isinstance(dsep_equations_to_run[0], dict):'
)

with open("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/api.py", "w") as f:
    f.write(new_content)
