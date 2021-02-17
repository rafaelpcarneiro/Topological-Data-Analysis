#  Extracting topological features from data 

I intend to keep here, at this repository, all of my work
concerning to the study of topological data analysis,
which will be called as TDA (to keep it short).

## So, what is TDA.

Basically, the idea here is to study data with a more
geometrically perspective. And this is done by relying
on the field of algebraic topology called homology.
Features such as connecivity, holes and etc will be analyzed
by a chain of homological groups, which will be indexed
by some filtration of topological spaces (that can be triangulated)
or simplicial complexes (which include graphs).

A good reference to study the subject can be obtained by the first three
chapters of the book:
* **Persistence Theory: From Quiver Representations to Data Analysis**,
  by the author Steve Oudot. 
This book as well other interesting texts can be found at the webpage
of [Steve Oudot](https://geometrica.saclay.inria.fr/team/Steve.Oudot/)

## A list of what I have done.
Here is a list of the programs or texts I have done so far.
Every work is listed down bellow and a brief notion of what
is done is given by some key words labelled in bold face
and enclosed with brackets.

1. **[Persistent Path Homology]**: Based on the paper, 
[**Persistent Path Homology of Directed Networks**](https://arxiv.org/abs/1701.00565) from
the authors Samir Chowdhury and Facundo Mémoli, I will implement
in python their algorithm responsible to calculate Persistent Path Homology 
from a network N=(X,A). This is done at the folder:
[persistent_path_homology](./persistent_path_homology)

2. 
