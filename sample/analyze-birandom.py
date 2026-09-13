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
# # Analysis for odd responses sqisign

# %%
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import math

# %%
p = 62

# Use the newest merged file for this p.
files = sorted(Path('.').glob(f'data/{p}-birandom-merged-*.csv.zip'))
if not files:
    raise FileNotFoundError(f'No merged CSV found for p={p}')
data_file = files[-1]
data_file

# %%
columns = (
    ['I1', 'I2', 'I3', 'I4']
    + ['OL2', 'OL3', 'OL4']
    + ['OR2', 'OR3', 'OR4']
    + [f'inner{i}{j}' for i in range(1, 5) for j in range(1, 5)]
)

df = pd.read_csv(data_file, compression='zip', comment='#', header=None, names=columns)
df

# %%
total_record = len(df)
total_record

# %% [markdown]
# We can easily check that it never happens that all 3 minima and their combinations are all even

# %%
print((df['inner11']*df['inner22']*df['inner12']*df['inner23']*df['inner13']*df['inner33']).sum()/total_record)

# %%
df[['I1', 'I2', 'I3', 'I4']].corr()

# %%
inner_cols = ['inner11', 'inner12', 'inner13', 'inner14',
   'inner21', 'inner22', 'inner23', 'inner24',
   'inner31', 'inner32', 'inner33', 'inner34',
   'inner41', 'inner42', 'inner43', 'inner44']
df["n_nonzero_inner"] = (df[inner_cols] != 0).sum(axis=1)
distribution = df["n_nonzero_inner"].value_counts().sort_index()
print(distribution/total_record)

# %%
inner_cols = ['inner11', 'inner12', 'inner13', 'inner22', 'inner23', 'inner33']
df["n_nonzero_inner"] = (df[inner_cols] != 0).sum(axis=1)
distribution = df["n_nonzero_inner"].value_counts().sort_index()
print(distribution/total_record)

# %%
df[inner_cols].corr()

# %%
(df['I3']).hist(bins=10000)

# %%
(df['I1']*df['I2']).hist(bins=1000)

# %%
df['I12'] = df['I1'] * df['I2']

# %%
df[['I12','OR2', 'OR3', 'OR4']].corr()

# %%
dfmin = df[['I1', 'I2', 'I3', 'I4']]
corr = dfmin.corr(method="spearman", numeric_only=True)
corr

# %%
df["I3-I1"] = df["I3"] - df["I1"]

# %% [markdown]
# ## 3rd minima estimates

# %%
d = (df["I3"]).to_numpy() 

n = len(d)
x = np.sort(d)
plt.scatter(np.arange(n)/n, x, s = 0.01)

# x values from 0 to max(d)
x_vals = np.linspace(0, d.max(), 500)

# y = number of samples bigger than x
counts = np.array([(d > x).sum() for x in x_vals]) / n

# optional: fit only the tail, e.g. fit_quantile=0.7
threshold = np.quantile(d, 1 - 0.99)





# %%

# x values from 0 to max(d)
lin_span = 1000
x_vals = np.linspace(0, d.max(), lin_span)

# y = number of samples bigger than x
counts = np.array([(d > x).sum() for x in x_vals]) / n

# optional: fit only the tail, e.g. fit_quantile=0.7
threshold = np.quantile(d, 1 - 0.5)

mask = (x_vals > 0) & (counts > 0) & (x_vals >= threshold) #& (x_vals < 20)


# %%
threshold_low = np.quantile(d, 0.2)
threshold_up = np.quantile(d, 1 - 0.00001)
print(threshold_low, threshold_up)
mask = (x_vals > 0) & (counts > 0) & (x_vals >= threshold_low) & (x_vals <= threshold_up)

print(f'Considering only {sum(mask)}/{len(mask)} values')
log_x = np.log(x_vals[mask])
log_y = np.log(counts[mask])

slope, intercept = np.polyfit(log_x, log_y, 1)

a_hat = -slope
C_hat = np.exp(intercept)

fit_y = np.zeros_like(x_vals, dtype=float)

mask_pos = x_vals > 0
fit_y[mask_pos] = C_hat * x_vals[mask_pos] ** (-a_hat)

mask_plot = mask_pos & (fit_y < 1)

plt.scatter(x_vals[mask], counts[mask],
            label="empirical probability: # {L > x}",
            s=0.1)

plt.plot(x_vals[mask_plot], fit_y[mask_plot],
         label=f"fit: {C_hat:.2f} / x^{a_hat:.2f}",
         color="red", alpha = 0.4)

plt.xlabel("x")
plt.ylabel("# samples with L > x / tot")
plt.axvline(x = threshold, color = 'green', linestyle="--", label=f'{threshold:f}')

plt.legend()
plt.show()

# %%
a = 4
xv = x_vals[mask] ** -a
lyv = counts[mask]

slope, intercept = np.polyfit(xv, lyv, 1)

fit_y = np.zeros_like(x_vals, dtype=float)

mask_pos = x_vals > 0
fit_y[mask_pos] = slope * x_vals[mask_pos] ** (-a) + intercept

mask_plot = mask_pos & (fit_y < 1)

plt.scatter(x_vals[mask], counts[mask],
            label="empirical probability: # {L > x}",
            s=0.1)

plt.plot(x_vals[mask_plot], fit_y[mask_plot],
         label=f"fit: {slope:.2f} / x^{a} + {intercept:.2f}",
         color="red", alpha = 0.4)

plt.xlabel("x")
plt.ylabel("# samples with L > x / tot")
plt.axvline(x = threshold, color = 'green', linestyle="--", label=f'{threshold:f}')

plt.legend()
plt.show()

# %%
from scipy.stats import norm
import math

d = np.sort(df["I3"].dropna().to_numpy())
n = len(d)

alpha0 = 4

rows = []
perc_span = np.linspace(0.0002, 0.005, 100)

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


# %% [markdown]
# ## 1st and 2nd minima product estimate

# %%
# d = (df["I1"]*df["I2"]).to_numpy() 
d = (df["I1"]).to_numpy() 
n = len(d)
x = np.sort(d)
plt.scatter(np.arange(n)/n, x, s = 0.01)

# %%
# x values from 0 to max(d)
x_vals = np.linspace(0, d.max(), 10000)

# y = number of samples bigger than x
counts = np.array([(d < x).sum() for x in x_vals]) / n

# optional: only consider vectors shorter than c * p 
# seems, that c = 0.01 start to give results
threshold = 0.3

# the reason for this sieve is that we expect that for shorter vectors the heuristic will better hold
# and we only care about the short case, corresponding to the long tail for lambda_3
mask = (x_vals > 0) & (counts > 0) & (x_vals <= threshold)

# TODO: understand the number theory behind this

# %%

log_x = np.log(x_vals[mask])
log_y = np.log(counts[mask])

slope, intercept = np.polyfit(log_x, log_y, 1)

a_hat = slope
C_hat = np.exp(intercept)

fit_y = np.zeros_like(x_vals, dtype=float)

mask_pos = x_vals > 0
fit_y[mask_pos] = C_hat * x_vals[mask_pos] ** (a_hat)

mask_plot = mask_pos & (fit_y < 1) & (x_vals <= threshold*1.1)

plt.scatter(x_vals[mask], counts[mask],
            label="empirical probability: # {L < x}",
            s=2, marker = 'd')

plt.plot(x_vals[mask_plot], fit_y[mask_plot],
         label=f"fit: {C_hat:.2f}  x^{a_hat:.2f}",
         color="red", alpha=0.9, linestyle="--", linewidth=0.5)

plt.xlabel("x")
plt.ylabel("# samples with L < x / tot")
plt.axvline(x = threshold, color = 'green', linestyle="--", label=f'{threshold:f}', alpha=0.2)

plt.legend()
plt.show()

print(f'Considering {sum([ 1 for m in mask_plot if m])} points')

# %%
a = 2
xv = x_vals[mask] ** a
lyv = counts[mask]

slope, intercept = np.polyfit(xv, lyv, 1)

fit_y = np.zeros_like(x_vals, dtype=float)

mask_pos = x_vals > 0
fit_y[mask_pos] = slope * x_vals[mask_pos] ** (a) + intercept

mask_plot = mask_pos & (fit_y < 1) & (x_vals <= threshold*1.1)

plt.scatter(x_vals[mask], counts[mask],
            label="empirical probability: # {L > x}",
            s=0.1)

plt.plot(x_vals[mask_plot], fit_y[mask_plot],
         label=f"fit: {slope:.2f}  x^{a} + {intercept:.2f}",
         color="red", alpha = 0.4)

plt.xlabel("x")
plt.ylabel("# samples with L > x / tot")
plt.axvline(x = threshold, color = 'green', linestyle="--", label=f'{threshold:f}')

plt.legend()
plt.show()

# %%
