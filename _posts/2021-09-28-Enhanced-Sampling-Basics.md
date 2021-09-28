---
title: Basics of enhanced sampling
date: 2021-09-28
permalink: /posts/2021/09-Enhanced-Sampling-Basics
excerpt_separator: <!--more-->
toc: true
tags:
  - Enhanced Sampling
---


# Generalized Ensemble Sampling

## Background

The potential energy surface of a complex system may have multiple local minimums which would cause serious problems, for example, when simulating at a low energy, the system could be trapped in one local minimum state, resulting in wrong simulation results.

**Generalized Ensemble Sampling** is one approach to overcome such problem. When performing molecular dynamics sampling, the space with high potential energy can hardly be reached in a simulation and if the simulation is done in a generalized ensemble, every state is reweighted by a non-Boltzmann probability weight factor to make the probability of sampling each state with a equal probability. 

Such method allows **random walk** in the whole phase space and lets the simulation to escape from any energy barrier and thus having a better sampling of the whole phase space.

---

## Theories

### Multicanonical algorithm and simulated tempering

Consider N particles with mass $m_k$ with coordinate vector and momentum vector $\bold{q} ,\bold{p}$. The Hamiltonian can be expressed as:

$$H(\bold{q},\bold{p})=K(\bold{p})+E(\bold{q})$$

The kinetic part can be written as:

$$K(\bold{p})=\sum_{k=1}^{N}\frac{\bold{p_k}^2}{2m_k}$$

In canonical ensemble at temperature T, each state in phase space can be expressed as $x \equiv (q,p)$. The weighted factor such state can be evaluated with a Boltzmann factor:

$$W_B(x;T)=e^{-\beta H(q,p)}$$

The momentum and position are decoupled in the expression of the the constitute of Hamiltonian, meaning they are not directly influencing each other. And because the statistical average of the kinetic energy can be considered as $\frac{3}{2}Nk_BT$, the above Boltzmann factor can be considered as depending on the potential energy only.

$$W_B(x;T)=W_B(E;T)=e^{-\beta E}$$

As a consequence, the probability distribution of potential energy can be given as proportional to the density of states times the Boltzmann weight.

$$P_B(E;T)\propto n(E)W_B(E;T)$$

When increasing E, the density of states increases rapidly while the Boltzmann factor decreases exponentially. As a result, the probability distribution function can be a **bell-like** function that has a maximum around average energy at Temperature T.

 

In the multicanonical ensemble, each state is weighted by a non-Boltzmann factor to make the probability distribution function a constant.

$$P_{mu}(E;T)\propto n(E)W_{mu}(E;T)\equiv constant$$

And the factor can be written as

$$W_{mu}(E)\equiv e^{-\beta_0 E_{mu}(E;T_0)}=\frac{1}{n(E)}$$

Here, an arbitrary reference temperature is chosen, and the so-called multicanonical potential is defined by:

$$E_{mu} (E,T_0)=k_BT_0ln(n(E))=T_0S(E)$$

As the density of states is usually not known in prior, the multicanonical weight has to be determined numerically by iterations of preliminary runs.

## Iteration workflow

### Metropolis Monte Carlo:

The transition probability is given by

$$w(x\to x')= min(exp(-\beta_0 \Delta E_{mu}),1)$$

where

$$\Delta E_{mu}=E_{mu}(E';T_0)-E_{mu}(E;T_0)$$

### Molecular dynamics

For molecular dynamics simulation, the simulation should be performed at specific temperature $T_0$. The equation of motion should be revised to:

$$\dot{ \bold{p}}_k=-\frac{\partial E_{mu}(E;T_0)}{\partial \bold{q_k}}=-\frac{\partial E_{mu}(E;T_0)}{\partial E}\bold{f_k}$$

As we are performing simulation in an energy fixed ensemble, the expansion work and other work which would cause system energy change can be omitted. The above EOM can be further written as:

$$\dot{\bold{p}}_k = \frac{T_0}{T(E)}\bold{f_k}$$

The weight factor is determined by short trail simulations iteratively. The initial set up of the parameters should pick a high initial temperature $T_0$, while the initial trail energy and weight factor can be set as Boltzmann style (normal simulation style):

$$\left\{\begin{array}{l}E_{\mathrm{mu}}^{(1)}\left(E ; T_{0}\right)=E \\W_{\mathrm{mu}}^{(1)}\left(E ; T_{0}\right)=W_{\mathrm{B}}\left(E ; T_{0}\right)=\exp \left(-\beta_{0} E\right)\end{array}\right.$$

With such trail simulation can we have a bell shaped energy distribution and the corresponding energy can be defined as $E_{max}$

$$E_{max}=<E>_{T_0}$$

For each iteration, assume we are taking $i$th iteration, the current weight denoted as $W_{mu}^{(i)}(E;T_0)=exp(-\beta_0E_{mu}^{(i)}(E;T_0))$. The probability distribution function can be obtained correspondingly. During the current simulation, update another value $E_{min}^{(i)}$, which is the lowest-energy throughout all the simulations up to date. 

The next iteration can be initialized with

$$E_{\mathrm{mu}}^{(i+1)}\left(E ; T_{0}\right)= \begin{cases}E, & \text { for } E \geq E_{\max }, \\ E_{\mathrm{mu}}^{(i)}\left(E ; T_{0}\right)+k_{\mathrm{B}} T_{0} \ln N^{(i)}(E)-c^{(i)}, & \text { for } E_{\min }^{(i)} \leq E<E_{\max } \\ \left.\frac{\partial E_{\mathrm{mu}}^{(i+1)}\left(E ; T_{0}\right)}{\partial E}\right|_{E=E_{\min }^{(i)}}\left(E-E_{\min }^{(i)}\right)+E_{\mathrm{mu}}^{(i+1)}\left(E_{\min }^{(i)} ; T_{0}\right), & \text { for } E<E_{\min }^{(l)}\end{cases}$$

Where $c^i$ is introduced to diminish the discontinuity. 

$$c^{(i)}=k_BT_0lnN^{(i)}(E_{max})$$

This procedure is repeated until the energy distribution becomes flat. When the convergence is reached, the $E_{min}$ is the global minimum.

Any physical quantity can be obtained its expectation value at any temperature:

$$\langle A\rangle_{T}=\frac{\sum_{E} A(E) N_{\mathrm{mu}}(E) e^{\beta_{0} E_{\mathrm{mu}}\left(E ; T_{0}\right)-\beta E}} {\sum_{E} N_{\mathrm{mu}}(E) e^{\beta_{0} E_{\mathrm{mu}}\left(E ; T_{0}\right)-\beta E}}$$

## Simulated Tempering

In ST method, the temperature is not picked a fix value. When the temperature is not fixed, the weight factor can be defined as:

$$W_{ST}(E;T)=e^{-\beta E+a(T)}$$

Here, the energy is not modified while an appropriate function $a(T)$  should be chosen so that the probability distribution function can be flat.

$$\begin{align} P_{ST}(T)&=\int dE\ n(E)W_{ST}(E;T)\\&=\int dE\ n(E)e^{-\beta E+a(T)}\\&=constant \end{align}$$

In simulated tempering, the temperature should be sampled uniformly, which means a random walk in temperature space should be achieved. And this in turn induces random walk in potential energy space.

In actual algorithm, the temperature is chosen to be ranged from a very low value to ensure the minimum is the global minimum, and to a sufficient high value to ensure the system would not be trapped in specific local minimum. For the use of programming, the temperature is discretized to M values.

Now the function $a(T)$ can be reduced to $a_m$ at specific temperature $T_m$. Once the parameters $a_m$s are determined, and initial configuration and temperature being chosen, a simulated tempering can be realized by:

1. A canonical MC or MD at specific temperature $T_m$ is carried out for a certain steps.
2.  The temperature is updated to the neighboring values with the configuration fixed. The transition probability of such temperature updating is given with a Metropolis algorithm

$$w(T_m \to T_{m\pm 1})=\{\begin{array}{l}1     \\exp(-\Delta)\end{array} $$

$$\Delta=(\beta_{m\pm1}-\beta_{m})E-(a_{m\pm1}-a_m)$$

Here, I think the simulation in NVT MD, the temperature of the system is not a stable one as only the temperature of the external heat reservoir is determined, and the temperature of the system is equilibrated with the external reservoir. So I think this method may not be stable?
