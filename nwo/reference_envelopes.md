
# envelopes

After identifying the most probable anchor set, anchor-set-constrained
(ASC) posterior decoding is used to define the bounds and scores of
each individual domain in the target sequence.

## definitions of envelope components

### residues in domain $d$

The marginal probability that residue $i$ is in domain $d$ is denoted
$P(i \in d)$. The anchor residue $i^0_d$ for domain $d$ has $P(i^0_d
\in d) = 1.0$ by construction. $P(i \in d)$ can be obtained by
marginalizing over M|I states in ASC decoding UP and DOWN sectors for
domain $d$, or more efficiently from cumulative ASC decoding
probabilities of having reached a begin B state in the UP sector or E
state in DOWN. Residues in a path between B and the anchor, and the
anchor and E, must be emitted by match or insert states in the
homologous region for domain $d$.

$$
   P(i \in d) = \left\{ \begin{array}{ll}
             \sum_{j=i^0_{d-1}}^{i-1} \rho_B(j)  & \mbox{: UP sector } i = i^0_{d-1}+1 \ldots i^0_d-1 \\
			 1                                   & \mbox{: anchor } i = i^0_d \\
			 1 - \sum_{j=i^0_d}^{i-1} \rho_E(j)  & \mbox{: DOWN sector } i = i^0_d+1 \ldots i^0_{d+1}-1 \\
			 0                                   & \mbox{: all other } i\\
             \end{array}\right.
$$			 


### homologous region
      
The "homologous region" $h^a \ldots h^b$ for domain $d$ is defined as the
region that is more likely homologous to the model than not.

$$
\begin{eqnarray*}
     h^a & = & \min \{ i : i \leq i^0_d \mbox{ and } P(i \in d) > 0.5 \} \\
     h^b & = & \max \{ i : i \geq i^0_d \mbox{ and } P(i \in d) > 0.5 \}
\end{eqnarray*}
$$

The threshold of 0.5 is chosen to prevent overlapping homologous
regions. The $h^a \ldots h^b$ coords are only used in defining the
alignable homologous region $i^a \ldots i^b$; they are not recorded in
the envelope data structure.

### envelope coords
 
The "envelope" $i^a \ldots i^b$ denotes the bounds of the inferred
alignment, defined as the most probable start/end positions within the
homologous region $h^a \ldots h^b$:

$$
\begin{eqnarray*}
      i^a & = & 1 + \mathrm{argmax}_{i=h^a-1}^{i^0_d-1} \rho_B(i) \\
      i^b & = &     \mathrm{argmax}_{i=i^0_d}^{h^b}     \rho_E(i)
\end{eqnarray*}
$$

This makes sure that glocal alignments end at plausible endpoints,
which probability mass alone (i.e. the homologous region) does not
assure.


### complete vs fragment domains

Whether domain $d$ is considered to be complete (glocal alignment) or
a fragment (local alignment) is determined from marginal probabilities
of using the G vs L state in the anchor-constrained path ensemble that
passes through domain $d$'s anchor, marginalized over the uncertain
start position:

$$
      p_L(d) = \sum_{i=i^0_{d-1}}^{i^0_d-1} \rho_L(i) \\
      p_G(d) = \sum_{i=i^0_{d-1}}^{i^0_d-1} \rho_G(i)
$$

Domain $d$ is marked as glocal if $p_G(d) \geq p_L(d)$, by setting its
`h4E_IS_GLOCAL` flag.

### outer envelope coords
 
 The "outer envelope" $o^a \ldots o^b$ includes the $i^a \ldots i^b$
 envelope, plus any additional residues with non-negligible
 probability mass. The outer envelope is used to calculate envelope
 scores (see below).  Outer envelope coords $o^a,o^b$ are defined the
 same way as the homologous region but with a lower threshold
 $\epsilon$:

$$
\begin{eqnarray*}
     o^a & = & \min \{ i : i \leq i^0(d) \mbox{ and } P(i \in d)  \geq \epsilon \} \\
     o^b & = & \max \{ i : i \geq i^0(d) \mbox{ and } P(i \in d)  \geq \epsilon \}
\end{eqnarray*}
$$

 Since epsilon < 0.5 (default 0.005), the outer envelope necessarily
 encompasses both the homologous region and the envelope:
 $o^a \leq h^a \leq i^a \leq i^b \leq h^b \leq o^b$.
        
### envelope score

The "envelope score" is the score of the ensemble of all paths through
this domain's anchor and no other, as if this were the only domain in
the entire sequence, constrained to start and end within the outer
envelope.
      
A domain is considered to be "distinct" if its outer envelope does not
overlap the outer envelope of an adjacent domain.  Then we know that:

$$
\begin{eqnarray*}
    \rho_{N|J}(o^a-1) & \geq & 1 - 2 \epsilon  \quad\mbox{and}\\ 
    \rho_{J|C}(o^b)   & \geq & 1 - 2 \epsilon 
\end{eqnarray*}
$$

on its left and right sides. When these conditions are met, we can
calculate the envelope score by a fast approximation using arithmetic
on existing ASC Forward matrix values, and we are guaranteed to obtain
a score within $\frac{4\epsilon}{\log 2}$ bits (default = 0.03 bits)
of the true envelope score. The `h4E_ENVSC_APPROX` flag is set for
domains where the approximation is used. Otherwise, if the domain is
not distinct, we recalculate the envelope score by ASC Forward on its
$o^a \ldots o^b$ outer envelope.
 


## off-by-one tricksiness notes

When we calculate $o^a,o^b$, we rearrange the definition of $o^b$ so
that both calculations are minima over i:

$$
\begin{eqnarray*}
    o^a_d & = & \min \{ i : i \leq i^0_d \quad\mbox{and}\quad P(i \in d) \geq \epsilon \} \\
    o^b_d & = & \min \{ i : i \geq i^0_d \quad\mbox{and}\quad P(i+1 \in	d)  < \epsilon \} 
\end{eqnarray*}
$$

Now the calculation can work in a single left to right pass of
increasing $i$, by identifying the first $i$ that satisfies the
criterion. We do the same for $h^a,h^b$.

This can get confusing, especially because of off-by-one issues. To think these
things through clearly, first think about what coords $i$ are even in
play. Domain $d$ must include its anchor $i^0_d$, and cannot include
the anchors of adjacent domains, so possible startpoints range from
$i^0_{d-1}+1$ to $i^0_d$, and possible endpoints range from $i^0_d$ to
$i^0_{d+1}-1$.
    
(Anchor sentinel conventions $i^0_0=0$ and $i^0_{D+1}=L+1$ make this
work even for first, last domains.)

Next, let's look at exactly how we have to do our sums over
$\rho_B(i)$ and $\rho_E(i)$.  The marginal probability that you end
exactly at residue $i$ is $\rho_E(i)$, but the probability that you
start exactly at $i$ is $\rho_B(i-1)$. Thus:

$$
  P(i \in d) = \left\{ \begin{array}{ll}
           \sum_{j=i^0_{d-1}}^{i-1} \rho_B(j)  & \mbox{:}\quad i = i^0_{d-1}+1 \ldots i^0_d \\       
     1.0 - \sum_{j=i^0_d}^{i-1}    \rho_E(j)  & \mbox{:}\quad i = i^0_d \ldots i^0_{d+1}-1
    \end{array} \right.
$$
           
That is, as you move left to right in $i$ on the left (starting) side
of the domain, the probability that $x_i$ is in the left side of the
domain is the sum of B's up to $i-1$, including B($i-1$) $\rightarrow$
M$_k(i)$ entry transitions at $x_i$. On the right (ending) side, the
probability that you're still in the domain starts at 1.0 at the
anchor $i^0_d$, and afterwards is 1.0 minus the sum of E up to i-1,
including M$_k(i-1) \rightarrow$ E exits on the previous row $i-1$.

Next, think about writing that down in pseudocode as two independent
calculations:

For $o^a$ start point:

```
pB = 0.0
for i = i0(d-1) to i0(d)-1   // note off by one!
   pB += B(i)                // pB is now P(i+1 \in d)
   if (pB >= epsilon)        // if i+1 has met threshold
     oa(d) = i+1             //   then i+1 is the first x_i satisfying it
     break
```

For $o^b$ end point:

```
pE = 0.0                      // pE starts with P(i0(d) \in d) = 1.0 at the anchor
for i = i0(d) to i0(d+1)-1:   // then for each possible endpoint i:
   pE += E(i)                 //   pE is now 1.0 - P(i+1 \in d)
   if (1.0 - pE < epsilon)    //   if i+1 has dropped below threshold
     ob(d) = i                //     then i was the last x_i satisfying it.
     break
```

Then the last thing to do is combine these two pieces in a single
calculation, over all domains and all i in a single pass, which is
what the implementation does.



## proof of envelope score approximation 

To prove that a domain is "distinct", so we can use a fast
approximation to obtain the envelope score, we want to identify that
the outer envelope bounds $o^a,o^b$ serve as "choke points" in the
N/C/J states, through which passes all but a negligible amount of path
mass.

Choke points don't have to exist, because path mass can flow into
$o^a+1 \ldots i^0$ from a previous domain, and out of $i^0 \ldots o^b-1$ to a next
domain; but if the domain is "distinct", these flows are
deemed negligible.

Let's make the notion of "negligible" more formal now.

From the domain number, we know whether the domain is flanked on the
left by N or J, and on the right by J or C; we call these X, Y
respectively. Let $q$ be the total amount of probability mass on paths
into the previous (resp. next) domain.

For $o^a$ as our left bound:

$$
\begin{array}{ll}
   P(o^a-1 \in d) < \epsilon               & \mbox{(1) by construction; that's how we chose } o^a. \\
   \rho_X(o^a-1) = 1 - P(o^a-1 \in d) - q  & \mbox{(2) mass flows back, accumulates in state X; eventually flows back to prev domain} \\
   1-\rho_X(o^a-1) < \epsilon + q          & \mbox{(3) substitute (1) in (2)}\\
   \rho_X(o^a-1) \geq 1 - \epsilon - q.    & \mbox{(4) rearrange (3)}\\
   \mbox{If $q \geq \epsilon$:}            & \mbox{(5)} \\
   \qquad\rho_X(o^a-1) \geq 1 - 2\epsilon  & \mbox{(6) substitute (5) in (4)}\\
\end{array}
$$
   
So, our definition of "distinct" at the left edge is that
if $\rho_X(o^a-1) \geq 1-2\epsilon$, we know that the amount of probability
mass in paths that enter to our right from a previous domain
(without bottlenecking through X(o^a-1)) is $\leq \epsilon$, and we know
that the amount of probability mass that's still in the homology
model (outside of the outer envelope) is $< \epsilon$. So we've lost
up to $2\epsilon$ of the probability mass in paths at the left edge.

Analogous analysis holds for right edge. 

Therefore, if the condition on X$(o^a-1)$ and Y$(o^b)$ hold, we know
that if we only look at paths that pass through X$(o^a-1)$ and
Y$(o^b)$, we are neglecting at most a relative mass of $4\epsilon$ in
other paths that use this domain's anchor.

For $\epsilon = 0.005$, we lose up to a fraction 0.02 of the total
mass of the ensemble we count, corresponding to a maximum LLR score
loss of ~0.03 bits; we consider this to be a negligible score
difference.
