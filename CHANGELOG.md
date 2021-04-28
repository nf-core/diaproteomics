# nf-core/diaproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v.1.2.4 - [27.04.21]

- Raise memory requirements for FDR estimation and output visualization

## v.1.2.3 - [26.04.21]

- Optional mzTab output
- Optional BioReplicate annotation for MSstats
- Template update to 1.13.3

## v.1.2.2 - [24.02.21]

- Fix AWS full test profile
- Enable XIC compression by sqMass file format
- Update to DIAlignR@2119587
- Update dependencies

## v.1.2.1 - [18.02.21]

- Fix RT alignment of libraries by swapping x and y coordinates of linear regression
- Option to select irts from 1st and 4th RT quantile only to avoid overfitting to the center of the RT distribution
- Update DIAlignR to 1.3.5@b0698a5
- Chromatogram indexing at low memory usage

## v.1.2.0 - [17.02.21]

- Skip DIA processing (if only library generation is needed)
- New DIAlignR parameters (optional parallelization)
- Update DIAlignR to 1.3.5@52eaf4e
- Multi-thread parameter to multiple steps
- Cache option
- Update tests
- Fixed TraML libray input
- Template update 1.12.1

## v.1.1.0 - [3.12.20]

- Template update
- DIAlignR 1.3.5 + parallel multicore execution
- Unify FDR scoring steps
- Add params

## v.1.0.0 - [13.11.20]

- 1.0.0 version
- fixed msstats input on peptide and protein level

## v1.0dev - [30.10.20]

- Change workflow step names
- Change column naming in input sheets
- Template update to 1.11
- Added full size test profile
- Added second test profile skipping library generation from DDA
- RawFile parsing for DDA input
- MSstats and general output plots
- Error message when using multiple spectral libraries but no merging
- Removed swath windows input

## v1.0dev - [15.09.20]

- Integrated sample sheet eg. experimental design
- Sample sheets can be supplied for dia and dda raw data, id and irt libraries
- more parameters
- RawFile parsing
- Pseudo iRT generation

## v1.0dev - [08.09.20]

- parameter documentation
- Template update 1.10.2

## v1.0dev - [24.07.20]

- integrated EasyPQP + testing

## v1.0dev - [15.07.20]

### `Added`

- Fresh start from Template 1.9
- Update of Container Dependencies
- Conceptualization for the use of EasyPQP for library generation and DIAlignR instead of TRIC

### `Deprecated`

- Initial workflow used TRIC alignment

## v1.0dev - [06.11.19]

- Initial release of nf-core/diaproteomics, created with the [nf-core](https://nf-co.re/) template.
