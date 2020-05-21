# Ficamacos

## Usage
Calculate option prices under Lévy market model using Filon's method, Carr-Madan's method and COS method.

The Carr-Madan's method using Carr-Madan's formula and fft follows from page 2-7 of http://homepages.ulb.ac.be/~cazizieh/sp_files/CarrMadan%201998.pdf

The Filon's method using Carr-Madan's formula follows from https://chasethedevil.github.io/lefloch_heston_filon.pdf page 2

The COS method using COS method given by Fang and Oosterlee https://link.springer.com/content/pdf/10.1007/s00211-009-0252-4.pdf


## Requirement
Numpy


## Installation

```
pip install Ficamacos
```

## Help
After installed the package, use following command to see the details of functions
```
from Ficamacos import model
from Ficamacos import method
help(model)
help(method)
```
The notation of parameters are the same as notation in the papers cited.
