/* C6-Tools Synthetic Biology Functions for Google Apps Script
    by J. Christopher Anderson with ChatGPT
    Copyright 2023 University of California, Berkeley
 */

/**
 * Translate a DNA sequence to its corresponding amino acid sequence.
 * 
 * @param {string} dna The DNA sequence to translate.
 * @return {string} The amino acid sequence encoded by the input DNA sequence.
 * 
 * @example
 * translate("ATGGAGAACTAG"); // returns "METEK*"
 */
// Create the genetic code object
  const geneticCode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
  };
function translate(dna) {
  // Ensure the input is a string
  if (typeof dna !== 'string') {
    throw new Error("Translate:  " + dna + " is not a string.");
  }
  
  // Ensure the input only contains valid DNA letters
  var validDNA = /^[ACGTacgt]+$/;
  if (!validDNA.test(dna)) {
    throw new Error("Input must only contain valid DNA letters (A, C, G, T).");
  }
  
  // Convert the input to uppercase
  dna = dna.toUpperCase();
  
  // Translate the DNA sequence to an amino acid sequence
  var aaSequence = "";
  for (var i = 0; i < dna.length; i += 3) {
        var codon = dna.substring(i, i + 3);
    var aa = geneticCode[codon];
    if (!aa) {
      throw new Error("Invalid codon: " + codon + "Input must be composed of valid triplets");
    }
    if(aa == '*') {
      continue;
    }
    aaSequence += aa;
  }
  
  return aaSequence;
}

/**
removeSites - A function that removes specified restriction enzyme sites from an open reading frame (ORF) DNA sequence
Usage:
removeSites(orf: string): string
@param {string} orf - A DNA open reading frame sequence
@return {string} A modified ORF DNA sequence where important restriction enzyme sites have been removed
Example:
removeSites("ATGTTCGGTCTCAACGGAGACCAGCAGGAATCTTAA");
Returns: "ATGTTCGGTCTGAACGGAGATCAGCAGGAATCTTAA"
*/
//codonUsageData maps amino acid codes to their most frequent codons in E. coli
//Codons with frequency < 0.1 are excluded
const codonUsageData = {
    F: ["TTT", "TTC"],
    S: ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    Y: ["TAT", "TAC"],
    C: ["TGT", "TGC"],
    L: ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    P: ["CCT", "CCC", "CCA", "CCG"],
    H: ["CAT", "CAC"],
    R: ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    Q: ["CAA", "CAG"],
    I: ["ATT", "ATC"],
    T: ["ACT", "ACC", "ACA", "ACG"],
    N: ["AAT", "AAC"],
    K: ["AAA", "AAG"],
    M: ["ATG"],
    W: ["TGG"],
    A: ["GCT", "GCC", "GCA", "GCG"],
    V: ["GTT", "GTC", "GTA", "GTG"],
    D: ["GAT", "GAC"],
    E: ["GAA", "GAG"],
    G: ["GGT", "GGC", "GGA", "GGG"],
};

//create array to hold forbidden sequences
	const forbiddenSequences = [];
	//loop over restriction enzymes
  for (const enzymeName in restrictionEnzymes) {
      const enzymeData = restrictionEnzymes[enzymeName];
      const recognitionSequence = enzymeData.recognitionSequence;
      forbiddenSequences.push(recognitionSequence);
      const rc = enzymeData.recognitionRC;
      if (recognitionSequence !== rc) {
          forbiddenSequences.push(rc);
      }
  }

function removeSites(orf) {
  //orf = "ATGTTCGGTCTCAACGGAGACCAGCAGGAATCTTAA";
	//validate input
  if (orf.length % 3 !== 0) {
      throw new Error("Invalid input sequence. Must be a multiple of 3.");
  }

  orf = orf.toUpperCase();
  if (!/^[ATGC]*$/.test(orf)) {
      throw new Error("Invalid input sequence. Must be composed of only DNA characters.");
  }

  //If there is a stop codon, preserve it, otherwise assume it is TAA
  var stopCodon = orf.substring(orf.length - 3);
  if (stopCodon === "TAA" || stopCodon === "TGA" || stopCodon === "TAG") {
      orf = orf.substring(0, orf.length - 3);
  } else {
    stopCodon = "TAA";
  }

	//split input into codon array
	const codonArray = orf.match(/.{1,3}/g);
	//translate orf to generate protein sequence
	const proteinSequence = translate(orf);

	outer: while (true) {
    var changeMade = false;
		//loop over forbidden sequences
		for (const forbiddenSequence of forbiddenSequences) {
			let searchStartIndex = 0;
			while (true) {
				//look for restriction site
				const enzymeIndex = orf.indexOf(forbiddenSequence, searchStartIndex);
				//if no site found, move on to next enzyme
				if (enzymeIndex === -1) {
					break;
				}

				//identify positions of codon array that overlap the restriction site
        var overlappingCodonIndices = [];
        for (var i = 0; i < Math.floor(forbiddenSequence.length/3); i++) {
            overlappingCodonIndices.push(i + Math.floor(enzymeIndex/3));
        }

        //randomly pick one of the overlapping codon indices
        const randomAAIndex = overlappingCodonIndices[Math.floor(Math.random() * overlappingCodonIndices.length)];
        const chosenCodon = codonArray[randomAAIndex];
        //get the amino acid code of the random codon from the translated ORF
        const aminoAcid = proteinSequence[randomAAIndex];
        //look up available codons for the amino acid
        const availableCodons = codonUsageData[aminoAcid];
        
        //randomly pick an earlier codon from the available codons
				let newCodon = availableCodons[Math.floor(Math.random() * availableCodons.length)];
				while (newCodon === chosenCodon) {
					newCodon = availableCodons[Math.floor(Math.random() * availableCodons.length)];
				}
        
				//replace the random codon in the codon array with the new codon
				codonArray[randomAAIndex] = newCodon;
        orf = orf.substring(0, (randomAAIndex * 3)) + newCodon + orf.substring((randomAAIndex * 3) + 3);
				searchStartIndex = enzymeIndex + 1;
        changeMade = true;
        continue outer;
			}
		}

    if(!changeMade) {
      break outer;
    }
	}  //end outer
	//rejoin the codon array to form the modified dna sequence
	return codonArray.join("") + stopCodon;
}

/**
 * Converts a peptide sequence of amino acids into a string of codons using the "one amino acid one codon" algorithm.
 *
 * The "one amino acid one codon" algorithm maps each amino acid in the input peptide sequence to a single codon.
 * This function uses the most frequently used codon for that amino acid in E. coli.
 *
 * @param {string} peptide - A string consisting of only amino acid letters and the asterisk character (*).
 * @return {string} A string of codons corresponding to the input peptide sequence.
 *
 * @throws {Error} If the input is not a string consisting of only amino acid letters and the asterisk character (*).
 *
 * @example
 * let peptide = "MFGLNGDQQES";
 * let codons = oneAAoneCodon(peptide);
 * console.log(codons);
 * // Output: "ATGTTTGGTTTAAATGGTGATCAACAAGAATCT"
 */
function oneAAoneCodon(peptide) {
    if (!/^[A-Z\*]+$/.test(peptide)) {
        throw new Error("Input must be a string consisting of only amino acid letters and the asterisk character (*)");
    }

    let codons = [];
    for (let i = 0; i < peptide.length; i++) {
        let aminoAcid = peptide[i];
        if(aminoAcid === '*') {
          codons.push('TAA');
        } else {
          let selectedCodon = codonUsageData[aminoAcid][0];
          codons.push(selectedCodon);
        }
    }
    return codons.join("");
}
