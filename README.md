# Rep-G-SixJSolver: A tool to Compute 6j Symbols of Rep(G)
### Background
Rep(G) refers to the category of complex representations of a finite group G. Since Rep(G) is a unitary fusion category, it is equipped with associators. This repository contains GAP files that can explicitly compute these associators. The entries of the associators are called 6j symbols (which is why this project is named as a _Rep(G) 6j Solver_)

**Data on Fusion Rings associated to Rep(G) (with/without multiplcities)**: The Grothedieck ring of Rep(G), denoted Gr(Rep(G)), is a commutative fusion ring consisting of representatives of isomorphism classes of irreducible representations of G as labels, equipped with + and x operators defined as the direct sum and tensor product of representations respectively. Along the way to compute 6j symbols, it is natural to first compute Gr(Rep(G)). A compilation of Grothendieck rings of the representation category of certain small finite groups is available.

**Relationship of 6j Symbols with Anyons and Planar Diagrams:** Since Rep(G) is a unitary fusion category, the Grothendieck ring of Rep(G) can be considered as an _anyon model_ (i.e. a commutative fusion ring which admits a categorification into a unitary fusion category). With access to explcit 6j symbols, we can construct a planar diagrammatic algebra representing anyon fusion/splitting and perform _F moves_ on these diagrams.  

