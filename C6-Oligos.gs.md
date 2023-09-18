 /**
 * @file C6-Oligos.gs
 * @author J. Christopher Anderson with ChatGPT
 * @copyright 2023 University of California, Berkeley
 * @license See the LICENSE file included in the repository
 * @version 1.0.0
 * @module C6-Oligos
 * @description
 * This script provides a collection of functions for designing oligos,
 * gblocks, or gene synthesis sequences and wetlab instructions for
 * manipulating them.
 *
 * @requires C6-Utils
 * @requires C6-Seq
 */

/**
 * Calculates a score for an annealing sequence based on various criteria.
 *
 * The score is calculated based on the following criteria:
 *  - Whether the first and last base of the annealing sequence are G or C
 *  - The G/C content of the annealing sequence (between 50% and 65%)
 *  - The base balance of the annealing sequence (no missing bases, no excesses of one particular base)
 *  - The randomness of the annealing sequence (no long stretches of a single base, no large regions of G/C rich sequence and all A/T in other regions)
 *
 * @param {string} inseq - The annealing sequence to score
 * @return {number} The score of the annealing sequence, between 0 and 1
 */
function scoreanneal(inseq) {
  anneal = resolveToSeq(inseq);
  let score = 0;
  const maxPossibleScore = 5;
  
  // Check if first and last base are G or C
  if (anneal[0] == "G" || anneal[0] == "C") {
    score++;
  }
  if (anneal[anneal.length - 1] == "G" || anneal[anneal.length - 1] == "C") {
    score++;
  }
  
  // Check G/C content
  const gcContent = gccontent(anneal);
  if (gcContent >= 0.5 && gcContent <= 0.65) {
    score++;
  }
  
  // Check base composition
  const baseBalance = basebalance(anneal);
  if (baseBalance > 0.75) {
    score++;
  }
  
  // Check for randomness
  const maxRepeat = maxrepeat(anneal);
  if (maxRepeat <= 3) {
    score++;
  }
  
  // Check length
  const lengthDiff = Math.abs(anneal.length - 20);
  score -= lengthDiff/2;
  
  return Math.max(0,score / maxPossibleScore);
}

/**
 * Returns the best annealing sequence for an input DNA sequence.
 *
 * The annealing sequence is a substring of the input sequence that meets the following criteria:
 *  - It is between 18 and 25 bases in length
 *  - If lock5 is true, the start of the annealing sequence must be the start of the input sequence
 *  - If lock3 is true, the end of the annealing sequence must be the end of the input sequence
 *  - The first and last base of the annealing sequence are ideally G or C
 *  - The overall G/C content of the annealing sequence is between 50% and 65%
 *  - The base composition of the annealing sequence is balanced (no missing bases, no excesses of one particular base)
 *  - The annealing sequence appears random (no long stretches of a single base, no large regions of G/C rich sequence and all A/T in other regions)
 *
 * If both lock5 and lock3 are true, the function throws an Error. It is not topologically possible.
 *
 * @param {string} inseq - The input DNA sequence
 * @param {boolean} lock5 - Whether the annealing sequence must start at the start of the input sequence
 * @param {boolean} lock3 - Whether the annealing sequence must end at the end of the input sequence
 * @return {string} The best annealing sequence that meets the specified criteria
 */
function findanneal(inseq, lock5, lock3) {
  inseq = resolveToSeq(inseq);

  const minLength = 18;
  const maxLength = 25;
  let bestAnneal = "N/A";
  let bestScore = -1;

  // Case: lock5=true, lock3=false
  if (lock5 && !lock3) {
    let startIndex = 0;
    for (let endIndex = minLength; endIndex <= maxLength; endIndex ++) {
      let anneal = inseq.substring(startIndex, endIndex);
      let score = scoreanneal(anneal);
      if(score > bestScore) {
        bestAnneal = anneal;
        bestScore = score;
      }
    }
    return bestAnneal;
  }

  // Case: lock5=false, lock3=true
  if (!lock5 && lock3) {
    let endIndex = inseq.length;
    for (let startIndex = endIndex - maxLength; startIndex < endIndex - minLength; startIndex++) {
      let anneal = inseq.substring(startIndex, endIndex);
      let score = scoreanneal(anneal);
      if(score > bestScore) {
        bestAnneal = anneal;
        bestScore = score;
      }
    }
    return bestAnneal;
  }

  if (!lock5 && !lock3) {
    let annealStart = 0;
    let annealEnd = inseq.length;
    
    for (; annealStart < annealEnd - minLength; annealStart++) {
      for (let i = annealStart + minLength; i < annealEnd; i++) {
        const anneal = inseq.substring(annealStart, i);
        const score = scoreanneal(anneal);
        
        if (score > bestScore) {
          bestAnneal = anneal;
          bestScore = score;
        }
      }
    }
    
    return bestAnneal;
  }

  throw new Error(`Cannot lock both ends of the template`);
}

/**
 * PCA function that generates oligos for a synthon using the polymerase chain assembly method.
 *
 * @param {string} synthon - The synthon DNA sequence.
 * @return {string} - A JSON array of oligos needed to build the synthon.
 */
function pca(synthon) {
  synthon = resolveToSeq(synthon);

  //Figure out the spacing
  var chunks = Math.round(synthon.length / 25);
  if(chunks % 2 != 0) {
    chunks++;
  }
  const chunksize = Math.round(synthon.length/chunks);

  // Initialize an array to store the annealing indices
  const annealingIndices = [];

  // Iterate through the synthon chunksize (about 25) bp at a time
  for (let i = chunksize; i < synthon.length - chunksize; i += chunksize) {
    // Get the 12 bp before and after the current site
    const seq = synthon.substr(i - 12, 24);

    // Get the best annealing site within the current window
    const anneal = findanneal(seq, false, false);

    // Store the start and end indices of 'anneal' on 'synthon'
    const startIndex = synthon.indexOf(anneal);
    const endIndex = startIndex + anneal.length;
    annealingIndices.push([startIndex, endIndex]);  
  }

  // Initialize an array to store the annealing indices
  const oligos = [];

  // Iterate through annealingIndices and construct oligos
  for (let i = 0; i < annealingIndices.length; i++) {
    //The forward first oligo
    if(i==0) {
      oligos.push(synthon.substring(0, annealingIndices[1][1]));
      continue;
    }

    //The last reverse oligo
    if(i==annealingIndices.length - 1) {
      var lastoligo = synthon.substring(annealingIndices[i][0]);
      oligos.push(revcomp(lastoligo));
      continue;
    }

    //For internal forward oligos
    if(i%2 == 0) {
      oligos.push(synthon.substring(annealingIndices[i][0], annealingIndices[i+1][1]));
      continue;
    }

    //For internal reverse oligos
    var rcoligo = synthon.substring(annealingIndices[i][0], annealingIndices[i+1][1]);
    oligos.push(revcomp(rcoligo));
  }

  // Return the oligos as a JSON array
  return oligos;
}

/**
 * LCA function that generates oligos for a synthon using the ligase chain assembly method.
 *
 * @param {string} synthon - The synthon DNA sequence.
 * @return {string} - A JSON array of oligos needed to build the synthon.
 */
function lca(synthon) {
  synthon = resolveToSeq(synthon);
    var seqLen = synthon.length;
    var oligos = [];

    // Find the largest chunk size that allows for roughly equal-sized chunks
    var mod = (seqLen - 25) % 50;
    var n = (seqLen - mod) / 50;
    if (mod > 25) {
        n++;
    } else if (mod < -25) {
        n--;
    }
    var chunkSize = Math.floor(seqLen / n);

    function generateOligos(s) {
        // Generate the n chunkSize oligos
        for (var i = 0; i < n; i++) {
            var oligo = s.substring(i*chunkSize, i*chunkSize + chunkSize);
            oligos.push(oligo);
        }
    }

    // Reverse complement the synthon sequence
    var synthonRevcomp = revcomp(synthon);

    // Generate oligos for the forward and reverse strands
    generateOligos(synthon);
    generateOligos(synthonRevcomp);

    // Return the oligos as a JSON array
    return oligos;
}

/**
 * Designs oligos or gene synthesis sequences for a BlgBrick part based on the input parameters.
 *
 * @param {string} sequence - the DNA sequence for the internal guts of the part.
 * @param {string} frgs - specifies whether a forward (F) or reverse (R) oligo is returned,
 *                      or a sequence for gene synthesis synthon (S) or a linear gblock (G).
 *
 * @returns {string} - the designed sequence.
 *
 * @example
 *
 * // Design a forward oligo for a BlgBrick part
 * const partSequence = 'atgcatgtaagtaattttacagctggattgctattacttgtaatagcatttggcggaacataa';
 * const forwardOligo = bglbrick(partSequence, 'F'); // returns 'ccataAGATCTATGCATGTAAGTAATTTTAC'
 *
 * // Design a reverse oligo for a BlgBrick part
 * const reverseOligo = bglbrick(partSequence, 'R'); // returns 'catcaCTCGAGttaGGATCCTTATGTTCCGCCAAATGCTA'
 *
 * // Design a gene synthesis sequence for a BlgBrick part
 * const geneSynthesisSeq = bglbrick(partSequence, 'G'); // returns 'AGATCTggataGAATTCatgAGATCTATGCATGTAAGTAATTTTACGGATCCtaaCTCGAG'
 */
function bglbrick(sequence, frgs) {
  rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);

  if ( rORf === 'F') {
    return "ccata" + "AGATCT" + findanneal(sequence, true, false);

  } else if (rORf === 'R') {
    return "catca" + "CTCGAGttaGGATCC" + revcomp(findanneal(sequence, false, true));

  } else if (rORf === 'S') {  
    return "GAATTCatgAGATCT" + sequence + "GGATCCtaaCTCGAG";

  } else if (rORf === 'G') {  
    return "ccataGAATTCatgAGATCT" + sequence + "GGATCCtaaCTCGAGtaacg";

  } else {
    throw new Error("Invalid value for 'frgs'. Please enter either 'F' or 'R' or 'Gblock' or 'Synthon'.");
  }
}

/**
 * Designs oligos or gene synthesis sequences for a BioBrick (RFC10).
 *
 * @param {string} sequence - the DNA sequence for the internal guts of the part.
 * @param {boolean} isCDS - specifies if the sequence is a CDS (true) or not (false).
 * @param {string} frgs - specifies whether a forward (F) or reverse (R) oligo is returned,
 *                      or a sequence for gene synthesis synthon (S) or a linear gblock (G).
 *
 * @returns {string} - the designed sequence.
 *
 * @example
 *
 * // Design a forward oligo for a BioBrick part with a CDS sequence
 * const partSequence = 'atgcatgtaagtaattttacagctggattgctattacttgtaatagcatttggcggaacataa';
 * const isCDS = true;
 * const forwardOligo = biobrick(partSequence, isCDS, 'F'); // returns 'gacttGAATTCgcggccgctTCTAGGGATAGAATTCATGAGATC'
 *
 * // Design a forward oligo for a BioBrick part with a non-CDS sequence
 * const partSequence = 'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac';
 * const isCDS = false;
 * const forwardOligo = biobrick(partSequence, isCDS, 'F'); // returns 'gacttGAATTCgcggccgctTCTAGAgTCCCTATCAGTGATAGAG'
 *
 * // Design a reverse oligo for a BioBrick part
 * const partSequence = 'atgcatgtaagtaattttacagctggattgctattacttgtaatagcatttggcggaacataa';
 * const isCDS = true;
 * const reverseOligo = biobrick(partSequence, isCDS, 'R'); // returns 'catcaACTAGTaTTATGTTCCGCCAAATGCTA'
 *
 * // Design a gene synthesis sequence for a BioBrick part with a CDS sequence
 * const geneSynthesisSeq = biobrick(partSequence, isCDS, 'G'); // returns 'GAATTCgcggccgctTCTAGatgcatgtaagtaattttacagctggattgctattacttgtaatagcatttggcggaacataatACTAGT'
 */
function biobrick(sequence, isCDS, frgs) {
  rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);

  if ( rORf === 'F') {
    if(isCDS) {
      return "gacttGAATTCgcggccgctTCTAG" + findanneal(sequence, true, false);
    } else {
      return "gacttGAATTCgcggccgctTCTAGAg" + findanneal(sequence, true, false);
    }

  } else if (rORf === 'R') {
    return "catca" + "ACTAGTa" + revcomp(findanneal(sequence, false, true));

  } else if (rORf === 'G') {  
    if(isCDS) {
      return "ccataGAATTCgcggccgctTCTAG" + sequence + "tACTAGTagcggccgCTGCAGcatcg";
    } else {
      return "ccataGAATTCgcggccgctTCTAG" + sequence + "tACTAGTagcggccgCTGCAGcatcg";
    }

  } else if (rORf === 'S') {  
    if(isCDS) {
      return "GAATTCgcggccgctTCTAG" + sequence + "tACTAGTagcggccgCTGCAG";
    } else {
      return "GAATTCgcggccgctTCTAGAg" + sequence + "tACTAGTagcggccgCTGCAG";
    }
    
  } else {
    throw new Error("Invalid value for 'frgs'. Please enter either 'F' or 'R' or 'Gblock' or 'Synthon'.");
  }
}

const stickyEnds = {
    UC: ['TACT', 'AAGC'],
    TP: ['GCTT', 'AGTA'],
    SP: ['AATG', 'ACCT'],
    U: ['TACT', 'CATT'],
    C: ['AGGT', 'AAGC'],
    T: ['GCTT', 'AGCG'],
    P: ['GGAG', 'AGTA']
};

/**
 * moclo - a function to design a forward or reverse oligo for MoClo cloning
 *
 * MoClo (Modular Cloning) is a DNA cloning method used in synthetic biology. 
 * It is based on the use of Type IIS restriction enzymes to generate standardized
 * 4 bp sticky ends based on part type. The types are P for promoter parts, U for 
 * 5' UTR parts (RBS), C for CDS/ORF parts, T for terminator parts, and SP for 
 * secretion tags.  UC parts are rbs.CDS joined parts.  TP parts are joined parts
 * of a terminator followed by a promoter.
 *
 * @param {string} sequence - the DNA sequence
 * @param {string} partType - one of {UC, TP, P, U, C, T, SP}
 * @param {string} frgs - specifies whether a forward (F) or reverse (R) oligo is returned,
 *                      or a sequence for gene synthesis synthon (S) or a linear gblock (G).
 *
 * @return {string} - the designed sequence
 *
 * @example
 * moclo("tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac", "P", forward);
 * // returns: "ccataGGTCTCaGGAGTCCCTATCAGTGATAGAG"
 */
function moclo(sequence, partType, frgs) {
  rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);


  if ( rORf === 'F') {
    var sticky = stickyEnds[partType][0];
    return "ccata" + "GGTCTCa" + sticky + findanneal(sequence, true, false);

  } else if (rORf === 'R') {
    var sticky = stickyEnds[partType][1];
    return "catca" + "GGTCTCt" + sticky + revcomp(findanneal(sequence, false, true));

  } else if (rORf === 'S') {  
    var sticky5 = stickyEnds[partType][0];
    var sticky3 = stickyEnds[partType][1];
    return "GGTCTCt" + sticky5 + sequence + revcomp(sticky3) + "aGAGACC";

  } else if (rORf === 'G') {  
    var sticky5 = stickyEnds[partType][0];
    var sticky3 = stickyEnds[partType][1];
    return "ccataGGTCTCt" + sticky5 + sequence + revcomp(sticky3) + "aGAGACCtaacg";

  } else {
    throw new Error("Invalid value for 'frgs'. Please enter either 'F' or 'R' or 'Gblock' or 'Synthon'.");
  }
}

/**
 * Designs oligos for PCR and subsequent homology-based assembly of two DNA sequences
 * into one molecule. It takes in two DNA sequences as input and returns the designed
 * oligos as output. The function is useful for joining overlapping DNA fragments or
 * inserting a DNA fragment into a larger construct using homology-based assembly methods
 * such as SOEing, Gibson assembly, or yeast homologous recombination.
 *
 * Usage:
 * var seq1 = "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac";
 * var seq2 = "atgcatgtaagtaattttacagctggattgctattacttgtaatagcatttggcggaacataa";
 * var oligos = genejoin(seq1, seq2, "F");
 * // returns: "GTGATAGAGATACTGAGCACATGCATGTAAGTAATTTTAC"
 *
 * @param {string} fivePrimeSeq The 5' DNA sequence to be joined.
 * @param {string} threePrimeSeq The second DNA sequence to be joined.
 * @param {string} ForR whether a forward (F) or reverse (R) oligo is returned
 * @return {string} The designed oligos for homology-based joining of the input sequences.
 */
function genejoin(fivePrimeSeq, threePrimeSeq, ForR) {
  fivePrimeSeq = resolveToSeq(fivePrimeSeq);
  threePrimeSeq = resolveToSeq(threePrimeSeq);

  var anneal5 = findanneal(fivePrimeSeq,false,true);
  var anneal3 = findanneal(threePrimeSeq,true,false);
  rORf = ForR[0].toUpperCase();

  var forOligo = anneal5 + anneal3;

  if ( rORf === 'F') {
    return forOligo;

  } else if (rORf === 'R') {
    return revcomp(forOligo);
  }

  throw new Error("Invalid value for 'ForR'. Please enter either 'Forward' or 'Reverse'.  You put in: " + ForR);
}

/**
 * Designs a ribosome binding site library of MoClo UC type
 *
 * @param {string} orf - the open reading frame of the CDS being controlled
 * @param {string} utr - the native or other known and function 5' UTR for the orf
 * @param {string} frg - whether a forward (F) or reverse (R) oligo is returned,
 * or a sequence for gene synthesis (G)
 *
 * @return {string} - the designed sequence
 */
function rbslib(orf, utr, frg)   {
  orf = resolveToSeq(orf);
  utr = resolveToSeq(utr);

  //Check that the orf is a multiple of 3 (codons)
  if (orf.length % 3 !== 0) {
    throw new Error("Length of orf must be a multiple of 3");
  }

  //If it lacks a stop codon, make it TAA
  if (!["TAA", "TGA", "TAG"].includes(orf.substring(orf.length-3))) {
      orf += "TAA";
  }

  //Abort if it isn't an orf
  if (!["ATG", "GTG", "CTG", "TTG"].includes(orf.substring(0, 3))) {
    throw new Error("Start of CDS not a start codon");
  }

  //Make it start with ATG
  if (!orf.startsWith("ATG")) {
    orf = "A" + orf.substring(1);
  }

  rORf = frg[0].toUpperCase();

  if (rORf === 'R') {
    return "catca" + "GGTCTCt" + "AAGC" + revcomp(findanneal(orf, false, true));
  }
  
  // Figure out the PWM library
  let rbs = utr.substring(utr.length - 7);
  rbs = "NVWGGRD" + rbs;
  rbs = utr.substring(utr.length - 17, utr.length - 14) + rbs;

  if ( rORf === 'F') {
    return "ccata" + "GGTCTCa" + "TACT" + rbs.toLowerCase() + findanneal(orf, true, false);
  }  

  if (rORf === 'G') {  
    return "ccataGGTCTCt" + "TACT" + rbs.toLowerCase() + orf + "GCTT" + "aGAGACCtgatg";
  }

  if (rORf === 'S') {  
    return "GGTCTCt" + "TACT" + rbs.toLowerCase() + orf + "GCTT" + "aGAGACC";
  }
}  

