
## H4_PROFILE - the HMMER4 dual-mode local/glocal probability model

### Boundary cases

In the state transitions t[0..M], t[0] and t[M] are special.
In the emission distributions e[0..M], e[0] is special.

t[0] stores G transitions in TM. The TI and TD distributions are
unused, but set to valid normalized probs anyway.

| t[0]| used for | set to  |
|-----|----------|---------|
| TMM | t(G->M1) | t(G->M1)|
| TMI | -        | 0       |
| TMD | t(G->D1) | t(G->D1)|
| TIM | -        | 1       |
| TII | -        | 0       |
| TID | -        | 0       |
| TDM | -        | 1       |
| TDI | -        | 0       |
| TDD | -        | 0       |
    
t[M] stores transitions to E. All of these probabilities are
constants, determined by the profile architecture.  The last M and D
state are forced to transit to E.

|t[M] | used for | set to  |
|-----|----------|---------|
| TMM | t(Mm->E) | 1       |
| TMI | -        | 0       |
| TMD | -        | 0       |
| TIM | -        | 1       |
| TII | -        | 0       |
| TID | -        | 0       |
| TDM | t(Dm->E) | 1       |
| TDI | -        | 0       |
| TDD | -        | 0       |

e[0] is unused. It is set with e[0][0] = 1, e[0][a>0] = 0.


### Differences relative to HMMER3 P7_HMM and P7_PROFILE

* Plan 9.

* No explicit insert emissions. They are implicitly equal to the
  background.

