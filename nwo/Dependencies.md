
## overview

| group        | depends on group(s)    |
|--------------|------------------------|
| base         |                        |
| models       | base                   |
| using models | base, models           |
| DP matrices  | base, models           |
| reference DP | base, models, matrices |


### base 

| .c file   | within-group dependencies |
|-----------|---------------------------|
| `logsum`  |                           |
| `simdvec` |                           |
| `general` | `logsum`                  |


### models 

Depends on _base_.

| .c file      | within-group dependencies |
|--------------|---------------------------|
| `h4_mode`    |                           |
| `h4_prior`   |                           |
| `h4_profile` |                           |
| `h4_counts`  | `h4_profile`              |
| `h4_hmmfile` | `h4_profile`              |
| `h4_path`    | `h4_mode`, `h4_profile`, `h4_counts`  |


### construction

Depends on _base_, _models_.

| .c file        | within-group dependencies    |
|----------------|------------------------------|
| `emit`         |                              |
| `modelsample`  |                              |
| `modelstats`   |                              |
| `parameterize` |                              | 
| `standardize`  |                              |
| `vectorize`    |                              |
| `zigar`        |                              |



### advanced construction

Depends on _base_, _models_, _construction_.

| .c file        | within-group dependencies    |
|----------------|------------------------------|
| `eweight`      |                              |
| `build`        | `eweight`                    |



### DP matrices

Depends on _base_, _models_.

| .c file        | within-group dependencies    |
|----------------|------------------------------|
| `h4_filtermx`  |                              |
| `h4_refmx`     |                              |


### reference DP

Depends on _base_, _models_, _matrices_.

| .c file        | within-group dependencies    |
|----------------|------------------------------|
| `reference_dp` |                              |


### vector DP

Depends on _base_, _models_, _matrices_.

| .c file         | within-group dependencies    |
|-----------------|------------------------------|
| `ssvfilter*`    |                              |
| `vitfilter*`    |                              |


### hmmer itself

Depends on everything.

| .c file         | within-group dependencies    |
|-----------------|------------------------------|
| `cmd_build`     |                              |
| `hmmer`         | `cmd_*`                      |


