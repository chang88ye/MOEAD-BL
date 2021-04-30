# **MOEA/D-BL**


MOEA/D-BL is a bio-inspired algorithm for multi/many-objective optimisation through a bilevel decomposition framework. MOEA/D-BL is essentially a derivative of MOEA/D (Zhang & Li, 2007), but it differs in the way how subproblems are created. It borrows the idea of M2M (Gu et al, 2014) to decompose an MOP into a sequence of subMOPs, and subMOPs are further broken down into scalar subproblems. Importantly, MOEA/D-BL allows neighboring subMOPs to shar subproblems so that positive information from one subMOP can be propagated into another and finally to all the subMOPs.

![](MOEADBL.png)

