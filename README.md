<img src="https://upload.wikimedia.org/wikipedia/commons/7/75/Logo_UNLP.jpg" alt="Logo of the project" width="400" height="200" align="right">

# Efficient simulation of non-markovian systems with Max-Ent type states &middot; [![Build Status](https://img.shields.io/travis/npm/npm/latest.svg?style=flat-square)](https://github.com/licTomasPerez) [![npm](https://img.shields.io/npm/v/npm.svg?style=flat-square)](https://www.npmjs.com/package/npm) [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com) [![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://github.com/your/your-project/blob/master/LICENSE)


## Introduction and Theoretical background

In my licenciatura in Physics' (roughly Bachelor+Master in Physics equivalent) thesis , with my director Dr. J.M. Matera, we worked in both analytical and computational techniques for the study of non-markovian quantum systems using Max-Ent approximation for instantaneous states, the latter with open code and hand-made also.
</p>

One of the fundamental problems in Quantum Information Theory and Quantum Computing is the accurate representation of quantum composite systems, in particular: their's states and dynamics. Said composite systems are univoquely represented by mathematical objects, the density operators, which containt the information of all possible n-body correlations present. A notable exception are the Gaussian dynamics. In these dynamics, achievable for bosonic systems, the dynamics is closed over the set of Gaussian states. Said set of Gaussian states are parameterizable in terms of pairwise correlations, thus forming a 𝑛 or 𝑛^2-dimensional Riemannian differentiable manifold; with a metric given by Hilbert-Schmidt inner product of operators. This has motivated to search for generalizations of the Gaussian states to still bosonic systems.

In this work, a generalization based on the Max-Ent property of Gaussian states is proposed, in which we will consider families of states that maximize the entropy of the system (Max-Ent principle), under the restriction of fixing the mean values of a certain set of independent observables. Strategies to build approximations within this family, which represent arbitrary states, will then be discussed.

## Objective

As an application case, we will study the relative entropy between the states that results from the dynamics in the Dicke model with its corresponding estimates as MaxEnt states defined by their local values and two-body correlations.

* We'll compare rho(t) with its max-ent state associated with a base of observables,
* compare rho(t) with its projected state, using the sc corresponding to the initial state, associated to a base of observables,
* and compare rho(t) with its projected state, using the sc corresponding to the instantaneous state, associated to a base of observables.

## Installing / Getting started

The present code was written and tested in Jupyter Notebook using Python version 3.9.7 using the Anaconda Distribution extension. 

## Developing

### Built With

In this work, the numerical routines for the Max-Ent and Projection-based approximations were written using the following packages 

<ol>
  <li>Numpy, used for building linear algebra subroutines and performing scientific calculations,</li>
  <li>Scipy, used for numerical optimization,</li>
  <li>Pickle, used for serializing and de-serializing a Python object structure, </li>
  <li>Matplotlib, used for data and results visualization</li>
  <li>QuTip, used for calculating quantum systems' dynamics and evolutions</li>
</ol>

### Prerequisites

To  be completed

### Setting up Dev

To be completed

### Building

To be completed

### Deploying / Publishing

## Versioning

To be completed

## Configuration

To be completed

## Tests

To be Completed

## Style guide

To be completed

## Licensing

To be completed

















