
# Zigars

A zigar is a BASE64-encoded, binary-compressed CIGAR representation of
a single domain alignment.

A bit string is encoded in BASE64 text using URL-safe, filename-safe
"base64url" `[A-Z][a-z][0-9]-_` encoding (RFC4648).  No = padding
characters are required (a zigar need not be a multiple of 4
char). The bit string is padded out with 0's to a multiple of six
(BASE64 consumes six bits per character).

The bit string is composed of a series of elements where each element
consists of a binary prefix code for an M/D/I symbol followed by a
binary prefix encoded runlength. The first element can be M, D, or I,
requiring 1 or 2 bits; subsequent elements only require 1 bit to
encode a change to one of two possible new states.  (HMMER4 alignments
never start or end with I, but zigars are specified more generally in
case they're used elsewhere.)

Each runlength $r$, $r>0$, is encoded by subtracting one and encoding
the $r-1$ value in a generalized exponential Golomb prefix code:
exp-Golomb-4 for M runlengths, exp-Golomb-0 for D/I runlengths. These
codes were chosen empirically to minimize zigar length on typical
HMMER4 alignments. 

Element encoding is summarized in this table:

|[]()| 0> | M> | I> | D> |   r-1 coded as: | codelen(r=1)| codelen(r=100k) |
|----|---:|---:|---:|---:|----------------:|------------:|----------------:|
| M  |  1 |  - |  0 |  0 |  exp-Golomb-4   |      5      |     29          |
| I  | 00 |  0 |  - |  1 |  exp-Golomb-0   |      1      |     33          |
| D  | 01 |  1 |  1 |  - |  exp-Golomb-0   |      1      |     33          |

A prefix of all 0's is an incomplete (unterminated) code for an
element, so end-of-data can be detected in a BASE64-encoded binary string
padded with trailing 0's.

No runlength exceeds 100K, given HMMER4 maximum profile and target
sequence lengths of 100K, so the maximum size of an element is 35
bits.  The implementation takes advantage of this by using just a
single 64bit unsigned integer as working space during encoding.

As a rule of thumb, an encoded element (state and runlength) requires
about 7.5 bits per M, 4.4 bits per I, and 3.6 bits per D. A typical
alignment has about 5 elements, requiring about 32 bits (albeit with
large variation). BASE64 consumes 6 bits per text character, so a
typical (mean) zigar is 5-6 characters long.






