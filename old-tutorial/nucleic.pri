# nucleic.pri
#
# Example of a prior file for DNA/RNA models.
# The values in this file are the HMMER 2 default settings.

Dirichlet   	# Strategy (mixture Dirichlet)
Nucleic		# type of prior (Amino or Nucleic)

# Transitions
1                     # Single component
1.0                   #   with probability = 1.0
0.7939 0.0278 0.0135  # m->m, m->i, m->d alpha's
0.1551 0.1331         # i->m, i->i alpha's 
0.9002 0.5630         # d->m, d->d alpha's

# Match emissions
# The use of 1.0 for alpha's here makes a simple Laplace "plus-one" prior.
#
1	# single component
1.0     #   with probability = 1.0
1.0 1.0 1.0 1.0

# Insert emissions
#  
1                   # Single component
1.0                 #    with probability 1.0
1.0 1.0 1.0 1.0
