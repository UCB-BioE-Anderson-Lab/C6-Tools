# C6-Tools v1.0 Synthetic Biology Functions for Google Apps Script

**by J. Christopher Anderson with ChatGPT**

_Last modified: June 12, 2023_

## Description

C6-Tools is a collection of JavaScript pure functions for use in Google Sheets to automate the design of synthetic biology experiments.

## Organization

It is organized into several script files:

- **C6-Utils**: Functions for manipulating JSON objects in a Sheet and their interconversions with ranges of worksheet cells.
  
- **C6-Seq**: Sequence utility functions and sequence operations. This handles the representation of DNA as strings or polynucleotide objects. It also has sequence operations including reverse complement (revcomp) and the like.

- **C6-Oligos**: Design of oligos for common DNA fabrication methods. Includes functions to design oligos (or gene synthesis cassettes) for MoClo Golden Gate Assembly, BioBricks, BglBricks, PCA, and LCA. It also has a lower level method, findanneal, for choosing annealing sequences during design of oligos for non-standardized experiments.

- **C6-Sim**: Simulation of molecular biology reactions. Simulates the products of PCR, Golden Gate Assembly, Gibson Isothermal Assembly/Yeast recombination. It can also parse and simulate construction files.

- **C6-Gene**: Design of coding sequences for gene synthesis. It can translate a DNA sequence to protein. It can silently remove restriction sites from a CDS. It can also reverse translate a protein sequence to DNA.

## Dependencies

None. This is a JavaScript project with no dependencies.

## Usage

To use these functions, the easiest method is to make a copy of a Google Sheet containing the scripts:

[C6-Tools v1.0](https://docs.google.com/spreadsheets/d/18GhA2s-x9kX1ar5YRMghXjcOHNW63eaZiD0fTT4xC60/edit?usp=sharing)

Alternatively, you could open a Google Sheet and navigate to the Tools menu. Select "Extensions > Apps Script" and paste in the scripts that you need. Save the scripts and return to your sheet. You can then call the functions from within any cell.

The distributions contain example Sheets showing how to use the functions for various tasks. Additional Sheets are used for testing and more nuanced usage illustration of some functions. You do not need to retain these examples for the tools to function properly; the Sheets can be deleted. The header comments for each function also include more technical information about their respective APIs.

## Limitations

Due to limitations in the Google Apps Script environment, these functions may be slower than equivalent functions running natively on a local computer. The degree of testing is uneven across this project, so beware that there may be bugs and errors.

## Contact

For questions, feedback, or report errors, contact [jcanderson@berkeley.edu](mailto:jcanderson@berkeley.edu).

## Acknowledgments

C6 is an homage to the Clotho effort pioneered with Douglas Densmore at UC Berkeley [Clotho History](https://clothocad.sourceforge.net/wiki/index.php/ClothoClassicHistory). There are 3 distinct releases of this Java-based synthetic biology platform. Projects C4, C5, and C6 represent ongoing explorations of how to interface with synthetic biology tools, refine the scope of use, and iterate on algorithm improvement. C6 is the first in the series written in JavaScript. It was inspired by the recognition that a Google Sheet could be easily distributed and would be readily recognizable to new users. Additionally, the spreadsheet interface is useful for visualizing experiments involving many samples and is commonly used in research. C6 was also motivated by the release of ChatGPT which greatly facilitated the translation of extensive legacy code to the JavaScript required for Google Sheets Apps Script. The project also borrows many ideas from the ConstructionFileSimulator effort, for which many students have contributed. The C6-Gene algorithms are inspired by the GeneDesign tool and follows a similar logic.

## Copyright

Copyright 2023 University of California, Berkeley

## License

See the LICENSE file included in the repository
