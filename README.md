# MEM AEAD Source Code Package

## Warning
The cipher designs of this source code package are very new and still **lack extensive analysis**. Therefore, **do not use** them in your applications just now!



## About
This repository provides implementations for the AEAD modes **OPP**, **MRO** and **MRS** instantiated with a round-reduced [BLAKE2b](https://blake2.net/) permutation. All ciphers target a 256-bit security level.

The specification of the schemes together with many more information can be found at [https://eprint.iacr.org/2015/999](https://eprint.iacr.org/2015/999).


### Features
* **OPP:**
    - based on the tweakable Masked Even-Mansour (MEM) block cipher
    - requires nonce-uniqueness
    - 1-pass
    - fully parallelisable
    - constant-time
* **MRO:**
    - based on the tweakable Masked Even-Mansour (MEM) block cipher
    - fully misuse-resistant
    - 2-pass
    - fully parallelisable
    - constant-time
* **MRS:**
    - based on the Sponge construction
    - fully misuse-resistant
    - 2-pass
    - constant-time

### Performance with 4 BLAKE2b Rounds

Platform     | Impl. |   OPP |   MRO |  MRS
-------------|-------|-------|-------|------
Cortex-A8    | NEON  |  4.26 |  8.07 | 8.50
Sandy Bridge | AVX   |  1.24 |  2.41 | 2.55
Haswell      | AVX2  |  0.55 |  1.06 | 2.40


### Performance with 6 BLAKE2b Rounds

Platform     | Impl. |  OPP |   MRO |   MRS
-------------|-------|------|-------|------
Cortex-A8    | NEON  | 5.91 | 11.32 | 12.21
Sandy Bridge | AVX   | 1.91 |  3.58 |  3.87
Haswell      | AVX2  | 0.75 |  1.39 |  3.58


## Designers

* [Robert Granger](http://people.epfl.ch/242282)
* [Philipp Jovanovic](https://zerobyte.io/)
* [Bart Mennink](http://homes.esat.kuleuven.be/~bmennink/)
* [Samuel Neves](https://eden.dei.uc.pt/~sneves/)

Code written by Philipp Jovanovic and Samuel Neves.

##License
The MEM AEAD source code is released under the [CC0 license](https://creativecommons.org/publicdomain/zero/1.0/). The full license text is included in the file `LICENSE`.
