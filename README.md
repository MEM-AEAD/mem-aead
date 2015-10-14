# MEM AEAD Source Code Package

## Warning
The cipher designs of this source code package are very new and still **lack extensive analysis**. Therefore, **do not use** them in your applications just now!

## About
The specification of the MEM family of AEAD schemes can be found on [ePrint](TODO). This repository provides implementations for the following ciphers:

* **MRO:** fully misuse-resistant, high-performance, based on the tweakable Masked Even-Mansour (MEM) block cipher
* **MRS:** fully misuse-resistant, high-performance, based on the Sponge construction
* **OPP:** nonce-respecting, high-performance, based on the tweakable Masked Even-Mansour (MEM) block cipher

All schemes leverage a round-reduced [BLAKE2b](https://blake2.net/) permutation.

## Designers

* [Robert Granger](http://people.epfl.ch/242282)
* [Philipp Jovanovic](https://zerobyte.io/)
* [Bart Mennink](http://homes.esat.kuleuven.be/~bmennink/)
* [Samuel Neves](https://eden.dei.uc.pt/~sneves/)

Code written by Philipp Jovanovic and Samuel Neves.

##License
The MEM AEAD source code is released under the [CC0 license](https://creativecommons.org/publicdomain/zero/1.0/). The full license text is included in the file `LICENSE`.
