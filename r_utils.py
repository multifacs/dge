import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2. robjects.packages import importr


def run_r(path: str):
    # Read the R script file
    with open(path, 'r') as f:
        r_script = f.read()

    # Execute the entire script

    robjects.r(f"suppressMessages({{\n{r_script}}})")
