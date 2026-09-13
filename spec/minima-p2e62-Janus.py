# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Exploring minima of the Janus lattice
#
# This notebook explores the distribution of minima of homsets of supersingular elliptic curves with scalar level structure. To run it you will need the following python packages:

# %%
import pandas as pd, numpy as np, gzip
from matplotlib import pyplot as plt

# %% [markdown]
# The dataset was pre-generated. To generate a new one, install Sagemath and then run
#
# ```
# ./sample_mt.py <how_many> -p<char> -g<mod>
# ```
#
# This generates `<how_many>` random ideal classes of B_p,∞ with `p = <char>` and prints the 4 minima of the associated quadratic form, normalized by √p, and the associated Gram matrix modulo `<mod>`.

# %%
# parse the data into a Pandas dataframe
with gzip.open('p2e62_janus.csv.gz', 'rt') as f:
    # first line is a (useful) comment
    for i in range(1):
        print(next(f))
    df = pd.read_csv(f, names=[f'm{i}' for i in range(1,5)], dtype=float)

# %% [markdown]
# ## The minima
#
# Remember the first 3 minima are divided by $p^⅓ / 2$, while the 4th is divided by $p/4$

# %%
df.describe()

# %% [markdown]
# 1st minimum, familiar shape of distribution

# %%
df['m1'].hist(bins=1000, color='blue')

# %% [markdown]
# first 3 minima, in log scale

# %%
np.log2(df['m1']).hist(bins=1000, color='blue', alpha=0.5, log=True)
np.log2(df['m2']).hist(bins=1000, color='red', alpha=0.5, log=True)
np.log2(df['m3']).hist(bins=1000, color='green', alpha=0.5, log=True)

# %% [markdown]
# The 4th minimum is essentially p/4

# %%
np.log2(df['m4']).hist(bins=1000, color='blue', log=True)

# %% [markdown]
# Maybe not so surprising, there's a decent correlation between 1st and 3rd minimum

# %%
np.log2(df['m1']*df['m3']).hist(bins=1000, log=True)

# %%
df['i2'] = 1/df['m2']
df['i3'] = 1/df['m3']
df['prod'] = df['m1'] * df['m2'] * df['m3']
df.corr()

# %%
df.plot('m1', 'i3', kind='scatter', s=0.1)

# %% [markdown]
# Usual quantile analysis in geometric progressoin (ratio 2)

# %%
df.quantile([1-1/2**i for i in range(1,21)])

# %%
df.quantile([1-1/2**i for i in range(1,21)]).pct_change()

# %%
import pandas as pd, numpy as np, gzip
from matplotlib import pyplot as plt
from scipy.stats import norm
import math

# %%
d = np.sort(1/df['m1'].dropna().to_numpy()) 
n = len(d) 
alpha0 = 1.5
rows = [] 
perc_span = np.linspace(0.001, 0.1, 100) 

for perc in perc_span:
    k = math.ceil(perc*n)
    u = d[-k-1]
    tail = d[-k:]
    gamma_hat = np.mean(np.log(tail / u))
    alpha_hat = 1 / gamma_hat
    
    se = alpha_hat / math.sqrt(k)
    
    # Test H0: alpha = 4 using the Hill asymptotic normal approximation
    se0 = alpha0 / math.sqrt(k)
    z = (alpha_hat - alpha0) / se0
    p_value = 2 * (1 - norm.cdf(abs(z)))
    
    rows.append({
        "k": k,
        "perc": perc,
        "threshold": u,
        "alpha_hat": alpha_hat,
        "ci_low": alpha_hat - 1.96 * se,
        "ci_high": alpha_hat + 1.96 * se,
        "z_alpha_eq_alpha0": z,
        "p_value_alpha_eq_alpha0": p_value,
        "C_hat_alpha0": (k / n) * u**4,
    })

out = pd.DataFrame(rows)
out

# %%
title = 'Hill Estimator'
fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

ax = axes[0]
ax.plot(out["perc"], out["alpha_hat"], marker=".", linewidth=0.3)
ax.fill_between(out["perc"], out["ci_low"], out["ci_high"], alpha=0.2)
ax.axhline(alpha0, color="red", linestyle="--", linewidth=1, label=fr"$\alpha={alpha0}$")
ax.set_ylabel(r"$\hat\alpha$")
ax.set_title(title or "Hill tail diagnostics")
ax.legend()

ax = axes[1]
ax.plot(out["perc"], out["p_value_alpha_eq_alpha0"], marker=".", linewidth=0.3)
p = 0.05
ax.axhline(p, color="red", linestyle="--", linewidth=1, label=f"{ p = }")
ax.set_ylabel(fr"p-value for $\alpha={alpha0}$")
ax.set_ylim(-0.02, 1.02)
ax.legend()

ax = axes[2]
ax.plot(out["perc"], out["C_hat_alpha0"], marker=".", linewidth=0.5)
C_0 = 0.1
ax.axhline(C_0, color="red", linestyle="--", linewidth=1, label=f"{ C_0 = }")
ax.set_ylabel(fr"$\hat C$ assuming $\alpha={alpha0}$")
ax.set_xlabel("perc")
ax.legend()

plt.tight_layout()
plt.show()


# %%
