# nf-core/diaproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [30.10.20]

### `Added`

Added full size test profile
RawFile parsing for DDA input
MSstats and general output plots
Error message when using multiple spectral libraries but no merging

## v1.0dev - [15.09.20]

### `Added`

Integrated sample sheet eg. experimental design
Sample sheets can be supplied for dia and dda raw data, id and irt libraries
more parameters
RawFile parsing
Pseudo iRT generation

## v1.0dev - [08.09.20]

### `Added`

parameter documentation
Template update 1.10.2

## v1.0dev - [24.07.20]

### `Added`

integrated EasyPQP + testing

## v1.0dev - [15.07.20]

### `Added`

Fresh start from Template 1.9
Update of Container Dependencies
Conceptualization for the use of EasyPQP for library generation and DIAlignR instead of TRIC

### `Deprecated`

Initial workflow used TRIC alignment

## v1.0dev - [06.11.19]

Initial release of nf-core/diaproteomics, created with the [nf-core](https://nf-co.re/) template.
