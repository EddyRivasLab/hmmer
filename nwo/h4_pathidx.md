
# H4_PATHIDX - index of domain locations in a path

The `H4_PATH` structure is highly compressed and memory-efficient.  
Sometimes it's convenient to associate more information with a path,
so that we can more easily access individual domains in it.

An `H4_PATHIDX` is always paired with a corresponding `H4_PATH`.

```
typedef struct {
  int *ia, *ib;
  int *ka, *kb;
  int *za, *kb;
  int *is_glocal;
  int  D;
  int  L;
} H4_PATHIDX;
```

* Domains are indexed 1..D.
   - accordingly, all `[0]` values (`ia[0]`, etc) for `d=0` are
     unused, and initialized to 0's just so they're initialized.

   - in the edge case of a zero-length homology path, D=0.  (See
     `h4_path.md` for more info on zero-length homology paths. They
     only arise in model construction from MSAs.)

* `ia[d]` and `ib[d]` are the start|end coordinates of domain `d`
  in the target sequence: 1..L.
   - L can be determined from the path, so we do, and store it.
  
* `ka[d]` and `kb[d]` are start|end coords on the profile: 1..M.
   - M can't be determined from an `H4_PATH`. We only know that 
     $ M \geq \max_d k^b(d) $.
  
* `za[d]` and `zb[d]` are the start|end coords in the path: 0..Z-1.[^1]
   - `za[d]` is the L|G state.
   - `zb[d]` is the last Mk|Dk profile state in the domain.

* `is_glocal[d]` is TRUE|FALSE, for whether this domain starts with L|G.
  
  
[^1]: If you really want to get into the details: since `st[0]` is
      always N and `st[Z-1]` is always C in a path, `za[d]` and
     `zb[d]` will never be 0 or Z-1.[^2]
	  
[^2]: And since L is always[^3] followed by an Mk, and G is
      always followed by Mk|Dk, `za[d]` can't be Z-2 either; in fact
      its valid range is 1..Z-3.

[^3]: In a local zero-length homology edge case of a `N-L-(E)-C` path
      (see `h4_path.md`), L can be followed by C; but in that case
      we're setting D=0 anyway because a zero-length homology isn't a
      "domain", and we set no coords [^4]
	
[^4]: Yes, I am a big fan of David Foster Wallace, how did you guess?
