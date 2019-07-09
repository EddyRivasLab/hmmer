
## Definitions of symbols used

### Roman font, upper case (predominately model states)

| symbol | meaning |
|--------|---------|
| B      |  begin-domain state |
| C      |  C-terminal (or 3') nonhomology state |
| D      |  delete state. DG = on glocal path, DL = local.  |
| E      |  end-domain state |
| G      |  glocal path begin state  |
| I      |  insert state. IG = on glocal path, IL = local.  |
| J      |  inter-domain (``join'') nonhomology state |
| L      |  local path begin state |
| M      |  match state. MG = on glocal path, ML = local. |
| N      |  N-terminal (or 5') nonhomology state |
| R      |  null (random) model's single state |
| S      | start state |
| T      | terminal state |


### math font, lower case (array indices, model parameters)

| symbol      | meaning |
|-------------|---------|
| $a$         | index in alphabet $1..K$ |
| $b$         | index in alphabet $1..K$ |
| $d$         | index for domains $1..D$ |
| $e$         | emission probability $e(k,a)$ |
|  ..         | (also: exponentiation, $e^x$)  |
| $f$         | null model emission probability $f_a$ |
| $g$         | index for segments in a sparse mask |
| $i$         | index on sequence $1..L$ |
| .. $i^a(d)$ | start of domain $d$ on sequence |
| .. $i^0(d)$ | anchor point for domain $d$ in sequence |
| .. $i^b(d)$ | end of domain $d$ on sequence |
| $j$         | index on sequence $1..L$  |
| $k$         | index on model $1..M$ |
| .. $k^a(d)$ | start of domain $d$ on model |
| .. $k^0(d)$ | anchor point for domain $d$ in model |
| .. $k^b(d)$ | end of domain $d$ on model |
| $l$         | _unused. too easily confused with 1._ |
| $o$         | outer envelope coordinates for a domain |
| .. $o^a(d)$ | outer envelope start for domain $d$ |
| .. $o^b(d)$ | outer envelope end for domain $d$ |
| $p$         | _probability of something_ |
| $s$         | scores |
| .. $s^F$    | a Forward score |
| .. $s^V$    | a Viterbi score |
| .. $s^A$    | an ASC Forward score |
| $t$         | transition probability $t_k(y)$ |
| $x$         | sequence   $x_1..x_L$ |
| $y$         | index/code for a state |
| $z$         | index in a path $\pi_{1..Z}$ |


### math font, upper case (includes sizes of arrays)

| symbol          | meaning |
|-----------------|---------|
| $D$             | number of domains identified in a sequence |
| $H$             | homology hypothesis (profile and its parameters) |
| $K$             | alphabet size  (typically 4 or 20) |
| $L$             | sequence or alignment length |
| $M$             | profile length; number of consensus match states  |
| $N$             | number of sequences |
| $P$             | _probability of something_ |
| $R$             | null hypothesis (null model and its parameters) |
| $S$             | output scores (bits) ... idealized |
| .. $S^{\dagger}$  | ... null2-corrected |
| .. $S^{\ddagger}$ | ... null2-corrected and scaled to units of bits |
| $W$             | null2 hypothesis (null2 model and its parameters) |
| $Z$             | path length |


### greek symbols, lower case

| symbol             | meaning |
|--------------------|---------|
| $\mathbf{\alpha}$  | Forward dynamic programming matrix |
| $\mathbf{\beta}$   | Backward dynamic programming matrix |
| $\mathbf{\gamma}$  | Viterbi dynamic programming matrix |
| $\delta$           | Delta function for maximum gain estimation alignment |
| $\epsilon$         | Loss parameter for ``well defined'' envelope definition; default = 0.005 |
| $\mathbf{\eta}$    | Null (nonhomology) model and its parameters |
| $\mathbf{\theta}$  | Profile and all its parameters |
| $\lambda$          | Label on a residue |
| \mathbf{\pi}$      | State path (traceback; alignment) |
| $\mathbf{\rho}$    | Posterior decoding dynamic programming matrix |
| $\sigma$           | Emission log-odds scores in the profile |
| $\tau$             | Transition scores in the profile |
| $\omega$           | Log-odds prior ratio for null2 vs. null hypothesis |


### greek symbols, upper case

| symbol             | meaning |
|--------------------|---------|
| $\Delta$           |log-odds ratio of priors; offset in converting score to posterior|



## Nomenclature conventions in the C code

### naming of constants

| prefix | module        | 
|--------|---------------|
| `h4_`  |  h4_profile, h4_mode, simdvec |
| `h4C_` |  h4_checkptmx |
| `h4D_` |  h4_domain    |
| `h4F_` |  h4_filtermx  |
| `h4H_` |  h4_hit       |
| `h4P_` |  h4_path      |
| `h4R_` |  h4_refmx     |
| `h4S_` |  h4_sparsemx  |




