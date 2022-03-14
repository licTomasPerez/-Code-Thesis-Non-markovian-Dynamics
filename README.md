# Efficient-simulation-of-non-markovian-dynamics-using-Max-Ent-states
In my licenciatura (roughly Bachelor+Master in Physics equivalent), with my director J.M. Matera, we worked in both analytical and computational techniques for the study of non-markovian quantum systems using Max-Ent approximation for instantaneous states, the latter with open code and hand-made also.

One of the fundamental problems in Quantum Information Theory and Quantum Computing is the accurate representation of quantum composite systems, in particular: their's states and dynamics. Said composite systems are univoquely represented by mathematical objects, the density operators, which containt the information of all possible n-body correlations present. A notable exception are the Gaussian dynamics. In these dynamics, achievable for bosonic systems, the dynamics is closed over the set of Gaussian states. Said set of Gaussian states are parameterizable in terms of pairwise correlations, thus forming a  𝑛  or  𝑛2 -dimensional Riemannian differentiable manifold; with a metric given by Hilbert-Schmidt inner product of operators. This has motivated to search for generalizations of the Gaussian states to still bosonic systems.

In this work, a generalization based on the Max-Ent property of Gaussian states is proposed, in which we will consider families of states that maximize the entropy of the system (Max-Ent principle), under the restriction of fixing the mean values of a certain set of independent observables. Strategies to build approximations within this family, which represent arbitrary states, will then be discussed.

As an application case, we will study the relative entropy between the states that results from the dynamics in the Dicke model with its corresponding estimates as MaxEnt states defined by their local values and two-body correlations.

* We'll compare rho(t) with its max-ent state associated with a base of observables,
* compare rho(t) with its projected state, using the sc corresponding to the initial state, associated to a base of observables,
* and compare rho(t) with its projected state, using the sc corresponding to the instantaneous state, associated to a base of observables.
