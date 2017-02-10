
## simdvec : config and initialization of SIMD vector environment


### on memory alignment, SIMD vectors, and malloc()

Yes, the C99 standard says malloc() must return a pointer "suitably
aligned so that it may be aligned to a pointer of any type of object"
(C99 7.20.3). Yes, __m128 vectors are 16 bytes. Yes, you'd think
malloc() ought to be required return a pointer aligned on a 16-byte
boundary. But, no, malloc() doesn't have to do this; malloc() will
return improperly aligned memory on many systems.  Reason: vectors are
officially considered "out of scope" of the C99 language.

So vector memory has to be manually aligned. For 16-byte alignment,
the idiom looks like the following, where for any to-be-aligned
ptr, we first allocate a raw (unaligned) memory space that's 15
bytes larger than we need, then set the aligned ptr into that space;
and when we free, we free the raw memory ptr:

```
  raw_mem = malloc(n + 15);
  ptr     = (__m128 *) ( ((uintptr_t) raw_mem + 15) & (~0xf));
```
  
More generally, for vector width W (in bytes):

```
  raw_mem = malloc(n + W - 1);
  ptr     = (__m256 *) ( ((uintptr_t) raw_mem + W - 1) & (~(W - 1)));
```

Technically, the result of casting a non-NULL pointer to an integer
type is undefined (C99 6.3.2.3), but this idiom for manual memory
alignment is in widespread use, and seems generally safe.

See also: posix_memalign(), as an alternative.
 
