<img src="https://upload.wikimedia.org/wikipedia/commons/7/75/Logo_UNLP.jpg" alt="Logo of the project" align="right">

# Efficient simulation of non-markovian systems with Max-Ent type states &middot; [![Build Status](https://img.shields.io/travis/npm/npm/latest.svg?style=flat-square)](https://travis-ci.org/npm/npm) [![npm](https://img.shields.io/npm/v/npm.svg?style=flat-square)](https://www.npmjs.com/package/npm) [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com) [![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://github.com/your/your-project/blob/master/LICENSE)
> Additional information or tag line

In my licenciatura (roughly Bachelor+Master in Physics equivalent), with my director J.M. Matera, we worked in both analytical and computational techniques for the study of non-markovian quantum systems using Max-Ent approximation for instantaneous states, the latter with open code and hand-made also.

One of the fundamental problems in Quantum Information Theory and Quantum Computing is the accurate representation of quantum composite systems, in particular: their's states and dynamics. Said composite systems are univoquely represented by mathematical objects, the density operators, which containt the information of all possible n-body correlations present. A notable exception are the Gaussian dynamics. In these dynamics, achievable for bosonic systems, the dynamics is closed over the set of Gaussian states. Said set of Gaussian states are parameterizable in terms of pairwise correlations, thus forming a  𝑛  or  𝑛2 -dimensional Riemannian differentiable manifold; with a metric given by Hilbert-Schmidt inner product of operators. This has motivated to search for generalizations of the Gaussian states to still bosonic systems.

In this work, a generalization based on the Max-Ent property of Gaussian states is proposed, in which we will consider families of states that maximize the entropy of the system (Max-Ent principle), under the restriction of fixing the mean values of a certain set of independent observables. Strategies to build approximations within this family, which represent arbitrary states, will then be discussed.

As an application case, we will study the relative entropy between the states that results from the dynamics in the Dicke model with its corresponding estimates as MaxEnt states defined by their local values and two-body correlations.

* We'll compare rho(t) with its max-ent state associated with a base of observables,
* compare rho(t) with its projected state, using the sc corresponding to the initial state, associated to a base of observables,
* and compare rho(t) with its projected state, using the sc corresponding to the instantaneous state, associated to a base of observables.



## Installing / Getting started

A quick introduction of the minimal setup you need to get a hello world up &
running.

```shell
commands here
```

Here you should say what actually happens when you execute the code above.

## Developing

### Built With
List main libraries, frameworks used including versions (React, Angular etc...)

### Prerequisites
What is needed to set up the dev environment. For instance, global dependencies or any other tools. include download links.


### Setting up Dev

Here's a brief intro about what a developer must do in order to start developing
the project further:

```shell
git clone https://github.com/your/your-project.git
cd your-project/
packagemanager install
```

And state what happens step-by-step. If there is any virtual environment, local server or database feeder needed, explain here.

### Building

If your project needs some additional steps for the developer to build the
project after some code changes, state them here. for example:

```shell
./configure
make
make install
```

Here again you should state what actually happens when the code above gets
executed.

### Deploying / Publishing
give instructions on how to build and release a new version
In case there's some step you have to take that publishes this project to a
server, this is the right time to state it.

```shell
packagemanager deploy your-project -s server.com -u username -p password
```

And again you'd need to tell what the previous code actually does.

## Versioning

We can maybe use [SemVer](http://semver.org/) for versioning. For the versions available, see the [link to tags on this repository](/tags).


## Configuration

Here you should write what are all of the configurations a user can enter when using the project.

## Tests

Describe and show how to run the tests with code examples.
Explain what these tests test and why.

```shell
Give an example
```

## Style guide

Explain your code style and show how to check it.

## Api Reference

If the api is external, link to api documentation. If not describe your api including authentication methods as well as explaining all the endpoints with their required parameters.


## Database

Explaining what database (and version) has been used. Provide download links.
Documents your database design and schemas, relations etc... 

## Licensing

State what the license is and how to find the text version of the license.

















