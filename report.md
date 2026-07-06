TL;DR: we need to understand how fast to the convergence from [Conjecture 6 of section 4.2](https://eprint.iacr.org/2024/1453.pdf), if we can get a bound of around $2^{-64}$ for a reasonable degree bound $B$ we are happy.



# Report for odd-degree

Recall: we want to understand how to pass from the isogeny sampler from SQIsign-odd
(that is SQIsign with odd-degree isogenies) to the isogeny sampler from
SQIsign-nist2 (that is SQIsign with any-degree isogenies).

## Notation

Let's first fix some notation.
We have and isogeny $\phi : E \to E'$, with $E$ informally the commitment
isogeny and $E'$ the challenge one.

We can define:

- $B \geq \sqrt{p}$, the bound on the degree of the isogeny $\phi$.
- $P_s$ is the set of pairs $(E,E')$ (commitment and challenge) that show up in
   SQIsign-nist2. We now that any of them has the same probability.
- $P_o$ is the set of pairs $(E,E')$ (commitment and challenge) that show up in
   SQIsign-odd, that is $(E,E') \in P_s$ such that there exists an 
   isogeny $\phi : E \to E'$ of odd degree and $deg(\phi) \leq B$.

Given $E,E'$ we can define:

- $B_s(E,E')$ is the set of isogenies $\phi : E \to E'$ of degree at most $B$.
- $B_o(E,E')$ is the set of isogenies $\phi : E \to E'$ of odd degree at most $B$.

## Assumption

From our analysis we can assume that there exists a constant $\epsilon_0 > 0$ such that

$| \frac{|B_o(E,E')|}{|B_s(E,E')|} - \frac{3}{8} | \leq \epsilon_0$

## Distributions

We can define the first two distributions we play with.

- $D_s$ is the distribution of isogenies in SQIsign-nist2. That is, we pick
  $(E,E') \in P_s$ uniformly at random and then pick $\phi \in B_s(E,E')$
  uniformly at random.
  So we have 
  $D_s(\phi) = \frac{1}{|P_s|} \frac{1}{|B_s(E,E')|}$

- $D_o$ is the distribution of isogenies in SQIsign-odd. That is, we pick
  $(E,E') \in P_o$ uniformly at random and then pick $\phi \in B_o(E,E')$ uniformly at random.
  So we have
  $D_o(\phi) = \frac{1}{|P_o|} \frac{1}{|B_o(E,E')|}$

We can now try to simulate the distribution $D_o$ using the distribution $D_s$.
The idea is to sample from $D_s$ and then reject samples that do not correspond
to odd-degree isogenies.
This give the distribution $D_{so}$, that is the conditioned probability of
$D_s$ on the event that the isogeny is of odd degree:

$D_{so}(\phi) =  D_s(\phi) / D_s(\bigcup_{(E,E') \in P_o} B_o(E,E'))$

that is:

$\frac{D_s(\phi)}{\sum_{P_o} \sum_{B_o(E,E')} D_s(\phi)}$

The numerator is $D_s(\phi) = \frac{1}{|P_s|} \frac{1}{|B_s(E,E')|}$ and the
denominator is

$\sum_{(E,E') \in P_o} \sum_{\phi \in B_o(E,E')} D_s(\phi) =$ 

$\sum_{(E,E') \in P_o} \sum_{\phi \in B_o(E,E')} \frac{1}{|P_s|} \frac{1}{|B_s(E,E')|} =$ 

$\sum_{(E,E') \in P_o} |B_o(E,E')| \frac{1}{|P_s|} \frac{1}{|B_s(E,E')|} =$ 

$\frac{1}{|P_s|} \sum_{(E,E') \in P_o}  \frac{|B_o(E,E')|}{|B_s(E,E')|} = \frac{s_o}{|PP_s|}$

with $s_o = \sum_{(E,E') \in P_o}  \frac{|B_o(E,E')|}{|B_s(E,E')|}$.


## Renyi divergence

If we compute the Renyi divergence of order $\infty$ between $D_o$ and $D_{so}$ we get:

$R_\infty(D_o || D_{so}) = \max_{\phi} \frac{D_o(\phi)}{D_{so}(\phi)} = \max_{\phi} \frac{\frac{1}{|P_o|} \frac{1}{|B_o(E,E')|}}{\frac{1}{|P_s|} \frac{1}{|B_s(E,E')|} / \frac{s_o}{|P_s|}} = \max_{\phi} \frac{|P_s| |B_s(E,E')| s_o}{|P_o| |B_o(E,E')| |P_s|} = \frac{s_o}{|P_o|} \max_{\phi} \frac{|B_s(E,E')|}{|B_o(E,E')|} =$ 

Using the definition of $\epsilon_0$ we get the upper bound: 

$\leq \frac{s_o}{|P_o|} (3/8 - \epsilon_0)^{-1} \leq \frac{3/8 + \epsilon_0}{3/8 - \epsilon_0}$.

Under the assumption that $0 \leq \epsilon_0 < 3/16$ (we can do better) we get the bound $\leq 1 + \frac{32}{3} \epsilon_0$.


## Security reduction


Let's assume that we have an adversary $\mathcal{A}$ that can break SQIsign-odd
with probability $\epsilon$.
Let $D_o(E)$ be the probability that $\mathcal{A}$ outputs a valid signature getting $q$ signing queries from the distribution $D_o$.
Let $D_{so}(E)$ be the probability that $\mathcal{A}$ outputs a valid signature getting $q$ signing queries from the distribution $D_{so}$.

By the [probability preservation of the Renyi divergence (Lemma 2.9)](https://eprint.iacr.org/2015/483.pdf) we have that:

$D_{so}(E) \geq D_o(E) / R_\infty(D_o \| D_{so})^q \geq D_o(E) / (1 + \epsilon)^q$
for $\epsilon = \frac{32}{3} \epsilon_0$.

Now, if we assume that $\epsilon \leq q^{-1}$ we have that $(1 + \epsilon)^q \leq e^{q \epsilon} \leq e^{1} \leq e$, so $D_{so}(E) \geq D_o(E) / e$.

In other words, if we have an adversary $\mathcal{A}$ that can break SQIsign-odd
with non-negligible probability $\kappa$ we can build an adversary $\mathcal{B}$ that can
break SQIsign-nist2 where we sieve for odd responses with still non-negligible probability at least $\kappa / e$.
Potentially using the Renyi divergence of order $a \approx \lambda$ we can get a better bound, but I am not sure, see more in the [relative error lemma](https://eprint.iacr.org/2017/480.pdf).


## Updates


A few updates on the previous analysis: there is an issue since it assumes that
the ration $\frac{|B_o(E,E')|}{|B_s(E,E')|}$ is the same for all $(E,E') \in
P_o$, but this is not true. In fact, there are several cases:

1. $\lambda_4 > B$, where on average the ratio is $3/8$.
2. $\lambda_3 < B < \lambda_4$, and the average ratio is $1/2$ or $1/4$ depending on the case.
3. $\lambda_2 < B < \lambda_3$, where the average ratio is $1/4, 1/2, 3/4$ depending on the case.
4. $\lambda_1 < B < \lambda_2$, where the average ratio is $1$.

The argument can be fixed, but we need to control some quantities per each of the cases listed above. Let $T$ be one of the cases, we can define:

- $P_{o,T} = \{ (E,E') \in P_o : (E,E') \text{ is of type } T \}$, the set of
  pairs $(E,E')$ that show up in SQIsign-odd and are of type $T$.
- $p_T$ average ratio of $\frac{|B_o(E,E')|}{|B_s(E,E')|}$ for $(E,E') \in
  P_{o,T}$.

We need to bound:

- maximum distance $\max_{(E,E') \in P_{o,T}} | \frac{|B_o(E,E')|}{|B_s(E,E')|} -  p_T |$ for each type $T$, this need to be approximately $2^{-64}$.
- ratio $\frac{|P_{o,T}|}{|P_o|}$ for each type $T$, that is expected to be
  close to $1$ for the first type and close to $0$ for the other types,


In particular, if this last ratio is approximately $2^{-64}$ for the non-first
types the rejection strategy that rejects all samples that are not of the first
type (i.e. $\lambda_4 > B$) will give the same previous reduction.

We can improve the strategy by also not rejecting samples for types $T$ that have a ratio $p_T$ larger or equal to $3/8$, but this seems a mess to implement.
It is interesting that this does not depends on the overall ration of pairs for
which no odd-degree isogeny exists.

I tried to adapt the previous analysis also to the rejection strategy that does
not check the type, using finite order Renyi divergence, but this is really tricky: it seems that we need a much smaller bound on the non-first types, in the order of $2^{-128}$, but this seems very strange, since it means that they would NEVER happen anyway in any universe.


Anyway, I have no idea on how to bound the ratios, so I will go to sleep for now xD

