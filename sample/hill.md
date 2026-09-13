# Hill Estimator Notes

These notes describe how we are using the Hill estimator to test the empirical claim

$$ \mathbb P(\lambda_3 > x) \sim \frac{C}{x^4} $$

for the third minimum of the random quaternionic lattice. In the data this quantity is the `I3` column.
Note that the data is scaled by $\sqrt{p}$, so the estimate of $C$ is for the
scaled variable $\lambda_3 / \sqrt{p}$.
The estimate for $\alpha$ is instead independent of the scaling.

## Tail Model

The Hill estimator is designed to estimate the behaviour of distributions with heavy upper tails. The standard model is
for distribution with tails behaving like the Pareto distribution, funny word to
say:

$$ \mathbb P(X > x) = \bar F(x) \sim C x^{-\alpha} $$

for large $x$, with $\alpha > 0$. 
In the extreme value theory usually they write instead using $$ \gamma = \frac{1}{\alpha}, $$

The Hill estimator estimates $\gamma$; so in our notation

$$ \hat\alpha = \frac{1}{\hat\gamma}. $$

For our case, $X = \lambda_3$, numerically represented by `df["I3"]`, and the hypothesis is $\alpha = 4$.

**Remark**:
Because the quaternionic setting has finitely many ideal classes, the
distribution is ultimately bounded: there is some large $K$ with
$\mathbb P(\lambda_3 > K) = 0$. Therefore the statement cannot be a literal
asymptotic as $x \to \infty$. The correct empirical statement is that, over the
observable range and well below the finite cutoff, the upper tail behaves like
an effective power law

$$ \mathbb P(\lambda_3 > x) \sim C x^{-4}. $$

## Definition

Consider the order statistics of the sample
$$
X_{(1)} \le X_{(2)} \le \cdots \le X_{(n)}
$$

(funny word to say we use the sorted observations).
Choose a tail size $k$, meaning that only the largest $k$ observations are used. The threshold is

$$
u = X_{(n-k)}.
$$

The Hill estimate of $\gamma = 1/\alpha$ is

$$
\hat\gamma(k)
= \frac{1}{k} \sum_{i=1}^k
\log\left(\frac{X_{(n-i+1)}}{u}\right).
$$

Then

$$
\hat\alpha(k) = \frac{1}{\hat\gamma(k)}.
$$

If the upper tail is close to Pareto with exponent $\alpha$, then $\hat\alpha(k)$ should be stable near $\alpha$ for a reasonable range of intermediate values of $k$.

## Choice of k

The parameter $k$ is the main tuning choice.

- If $k$ is too small, the estimate uses very few points and has high variance.
- If $k$ is too large, the estimate includes non-tail observations and becomes biased.
- The usual diagnostic is a Hill plot: plot $\hat\alpha(k)$ against $k$ and look for a stable plateau.

From the theory (Theorem 3.2.2 of de Haan, Ferreira) we should choose $k$ as a function of the sample size $n$ so that as $n \to \infty$
- $k(n) \to \infty$
- $k(n)/n \to 0$ .

In practice, we have a fixed sample size $n$, so we look for a range of $k$ that is not too small and not too large.



For our data, good diagnostics were obtained for values such that $k/n$ ranges between about $0.01 \%$ and $1 \%$ of the sample size (`n > 1e6`).
For larger k they started to behave differently, likely because it includes too much of the non-tail part of the distribution.

## Confidence Intervals

Under standard regular variation assumptions, the Hill estimator has an asymptotic normal approximation. For $\hat\alpha$, a rough standard error is

$$ \operatorname{se}(\hat\alpha) \approx \frac{\hat\alpha}{\sqrt{k}}. $$

Thus an approximate 95% confidence interval is

$$ \hat\alpha \pm 1.96 \frac{\hat\alpha}{\sqrt{k}}. $$

This interval should be treated as a diagnostic, not as a final rigorous theorem for the finite-support distribution. It is useful because the finite cutoff is far outside the sampled range.

## Testing $\alpha = 4$

For the null hypothesis

$$
H_0: \alpha = 4,
$$

we use

$$ z = \frac{\hat\alpha - 4}{4 / \sqrt{k}}. $$

The two-sided p-value is

$$
p = 2 \left(1 - \Phi(|z|)\right).
$$

This is again an asymptotic diagnostic. In our data, the p-values for $k$ between about $500$ and $20000$ did not reject $\alpha = 4$, while very large $k$ was less reliable.

## Estimating C

If we fix $\alpha = 4$, then at the threshold $u = X_{(n-k)}$ we have

$$
\mathbb P(X > u) \approx \frac{k}{n}
$$

and the model says

$$
\mathbb P(X > u) \approx \frac{C}{u^4}.
$$

Therefore

$$
\hat C(k) = \frac{k}{n} u^4.
$$

If the model $\mathbb P(X > x) \sim C/x^4$ is valid over the observed tail range, then $\hat C(k)$ should also be reasonably stable as a function of $k$.

## More todo

- [ ] apply this to different primes
- [ ] since it is normal we can do some bootstrapping of the data and see if the Hill estimator is stable under resampling


## References

Some references for the Hill estimator and extreme value theory, to be expanded.

1. Bruce M. Hill, "A Simple General Approach to Inference About the Tail of a Distribution", The Annals of Statistics 3(5), 1163-1174, 1975. DOI: 10.1214/AOS/1176343247.  
   https://doi.org/10.1214/AOS/1176343247

2. Laurens de Haan and Ana Ferreira, *Extreme Value Theory: An Introduction*, Springer, 2006. DOI: 10.1007/0-387-34471-3.  
   https://link.springer.com/book/10.1007/0-387-34471-3

3. `ExtremeRisks` R package documentation for `HTailIndex`.  
   https://search.r-project.org/CRAN/refmans/ExtremeRisks/html/HTailIndex.html
