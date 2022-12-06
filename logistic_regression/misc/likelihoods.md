Likelihoods
================
2022-11-18

**Likelihood (non-D&C)**

$$
L_i(\boldsymbol{\beta}) = \Big(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{y_i}\Big(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{1-y_i}
$$

**Multivariate normal prior on **β****

$$
\pi({\boldsymbol{\beta}}) \propto \exp\Big\\{ -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\Big\\}
$$

**Posterior (non-D&C)**

$$
\pi({\boldsymbol{\theta}}\|y_i) \propto \Big\[\Big(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{ y_i}\Big(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}\Big)^{1- y_i} \Big\] \cdot \Big\[\exp\Big\\{ -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\Big\\}  \Big\]
$$

*Log-posterior (non-D&C)*

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|y_i) ) &\propto y_i \log(\frac{\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}})+(1-y_i) \log(\frac{1}{1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\}}) -\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto  y_i\log(\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})- y_i\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})+ y_i \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
&\propto y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} -\log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})-\frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu})\\\\
\end{aligned}
$$

Therefore,

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|\mathrm{\bf{y}}) ) 
= \sum\_{i=1}^n (y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} - \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})) - \frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu}) + C\\\\
\end{aligned}
$$

for subset j

$$
\begin{aligned}
\log(\pi({\boldsymbol{\theta}}\|\mathrm{y}\_j) ) 
= \Big(\frac{n}{m_j}\Big) \sum\_{i=1}^{m_j} (y_i\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}} - \log(1+\exp\\{\mathrm{\bf{x}\_i^\top\boldsymbol{\beta}}\\})) - \frac{1}{2} (\boldsymbol{\beta}- \boldsymbol{\mu})^\top\Sigma^{-1}(\boldsymbol{\beta}- \boldsymbol{\mu}) + C\\\\
\end{aligned}
$$
