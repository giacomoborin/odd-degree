Note for comparison: the correct value for log of challenge degree in the even case is $\log p / 2 + \lambda$.
This because we need, by the soundness argument of SQIsign2D-west to have $e_{resp} + e_{ch} \leq e$, with $e$ the logarithm of the challenge degree.
For odd case is, given a failure $\leq 1/\epsilon$ the num of $(2,2)$ isogenies is $\log p / 2 + \log \epsilon / 4$.
For the case of:

| case | 2-iso | (2,2)-iso |
|------|-------|----------|
| even | $\log p / 2 + \lambda$ | $\log p / 2 $ |
| odd-$\log \epsilon = {\lambda}$ | $ \lambda$ | $\log p / 2 + \lambda/4$ |
| odd-$\log \epsilon = {64}$ | $ \lambda$ | $\log p / 2 + 16$ |
| odd-$\log \epsilon = {128}$ | $ \lambda$ | $\log p / 2 + 32$ |
| torsion | - | $\log p / 2 + \lambda/2 + 32$ |


The amortized difference (assuming one $(2,2)$-isogeny is equivalent to four $2$-isogenies) between the odd and even cases is $\log p /8 - \lambda/4$ for the first two cases.
We expect the odd degree variant to be roughly more efficient when $2 \lambda \leq  \log p$, that will be the case for the new parameters we will have to consider.

