# nf-core/diaproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v.1.2.3 - [14.04.21]

### `Added`

Optional mzTab output

## v.1.2.3dev - [25.03.21]

### `Added`

Template update to 1.13.3

## v.1.2.2 - [23.02.21]

### `Added`

Fix AWS full test profile
Enable XIC compression by sqMass file format
Update to DIAlignR@2119587
Update dependencies

## v.1.2.1 - [18.02.21]

### `Added`

Fix RT alignment of libraries by swapping x and y coordinates of linear regression
Option to select irts from 1st and 4th RT quantile only to avoid overfitting to the center of the RT distribution
Update DIAlignR to 1.3.5@b0698a5
Chromatogram indexing at low memory usage

## v.1.2.0 - [04.12.20]

### `Added`

Skip DIA processing (if only library generation is needed)
New DIAlignR parameters (optional parallelization)
Update DIAlignR to 1.3.5@52eaf4e
Multi-thread parameter to multiple steps
Cache option
Update tests
Fixed TraML libray input
Template update 1.12.1

## v.1.1.0 - [27.11.20]

### `Added`

Template update
DIAlignR 1.3.5 + parallel multicore execution
Unify FDR scoring steps
Add params

## v.1.0.0 - [11.11.20]

### `Added`

1.0.0 version
fixed msstats input on peptide and protein level

## v1.0dev - [30.10.20]

### `Added`

Change workflow step names
Change column naming in input sheets
Template update to 1.11
Added full size test profile
Added second test profile skipping library generation from DDA
RawFile parsing for DDA input
MSstats and general output plots
Error message when using multiple spectral libraries but no merging
Removed swath windows input

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
