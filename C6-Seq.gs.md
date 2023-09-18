/**
 * @file C6-Seq.gs
 * @author J. Christopher Anderson with ChatGPT
 * @copyright 2023 University of California, Berkeley
 * @license See the LICENSE file included in the repository
 * @version 1.0.0
 * @module C6-Seq
 * @description
 * This script provides a collection of functions for working with DNA sequences
 * and polynucleotide objects within Google Sheets. It offers various utilities
 * for sequence manipulation, validation, and analysis.
 *
 * @requires C6-Utils
 */

/**
 * @typedef {string} DNASequence
 * A string representing a DNA sequence with the following rules:
 * - Contains only the characters 'A', 'T', 'C', and 'G' (case-insensitive), or degeneracy codes
 * - No spaces or other characters are allowed
 */

/**
 * Validate and resolve a biopolymer sequence expressed as a string to a string of valid biopolymer letters
 * 
 * @param {string} sequence The biopolymer sequence to validate and resolve.
 * @return {string} The resolved and validated biopolymer sequence containing only valid biopolymer letters.
 * 
 * @example
 * cleanup("1 ATGGAGAACTAG GGTCTC"); // returns "ATGGAGAACTAGGGTCTC"
 * cleanup("AUGGAGAAACUAG\nGGUCUC"); // returns "AUGGAGAAACUAGGGUCUC"
 * cleanup("MVKHLIVTGLMVAL\nGLCSC"); // returns "MVKHLIVTGLMVALGLCSC"
 * @customfunction
 */
function cleanup(sequence) {
  // Ensure the input is a string
  if (typeof sequence !== 'string') {
    try {
      sequence = sequence.toString();
    } catch(err) {
      throw new Error("Input must be a string. It's a " + typeof sequence);
    }
  }
  
  // Remove any spaces, numbers, and line returns
  sequence = sequence.replace(/[\s\d\r]/g, "");
  
  // Ensure the input only contains valid biopolymer letters
  var validBiopolymer = /^[ACGTRYSWKMBDHVNUacgtryswkmbdhvnu*]+$/;
  if (!validBiopolymer.test(sequence)) {
    throw new Error("Input must only contain valid biopolymer letters (DNA: A, C, G, T, or degeneracy codes, RNA: A, C, G, U, or degeneracy codes, Protein: 20 standard amino acids and their degeneracy codes)");
  }
  
  // Convert the input to uppercase
  sequence = sequence.toUpperCase();
  
  return sequence;
}

/**
* This function takes in a string value, `seq`, and first converts it to a string if it is not. It then 
* matches the regular expression pattern defined by the constant _regexDNA.
* If the input matches the pattern, the input is returned.  Otherwise it throws an error.
* 
* @param {String} seq - The input containing a DNA sequence string.
* @return {String} - Returns the input sequence as a string
* 
* @example
* var sequence = resolveToSeq("ATCG");
* Logger.log(sequence); // outputs "ATCG"
* var sequence = resolveToSeq("name1");
* Logger.log(sequence); // outputs the column B value of the matching row
* @customfunction
*/
const _regexDNA = /^[ACTGactgMRWSYKVHDBNXmrwsykvhdbnx-]+$/;
function resolveToSeq(seq) {
  seq = seq.toString();
  if (_regexDNA.test(seq)) {
    return seq.toUpperCase();
  }

  throw new Error("Unrecognizable as sequence: " + seq);
}

/**
 * Represents a polynucleotide (DNA or RNA) molecule.
 *
 * @class
 * @param {string} sequence - The sequence of the polynucleotide.
 * @param {string} ext5 - The 5' extension of the current coding strand.
 * @param {string} ext3 - The 3' extension of the current coding strand.
 * @param {boolean} isDoubleStranded - Whether the polynucleotide is double stranded.
 * @param {boolean} isRNA - Whether the polynucleotide is RNA.
 * @param {boolean} isCircular - Whether the polynucleotide is circular.
 * @param {string} mod_ext5 - The 5' end modification.
 * @param {string} mod_ext3 - The 3' end modification.
 *
 * @example
 * var polynucleotide = new Polynucleotide("AGCTAGCT", "GATC", "CTAG", true, false, false, "Mod5", "Mod3");
 * console.log(polynucleotide.sequence); // Output: "AGCTAGCT"
 * console.log(polynucleotide.ext5); // Output: "GATC"
 * console.log(polynucleotide.isDoubleStranded); // Output: true
 */
class Polynucleotide {
  constructor(sequence, ext5 = null, ext3 = null, isDoubleStranded, isRNA, isCircular, mod_ext5, mod_ext3) {
    this.sequence = sequence.toUpperCase();
    this.ext5 = ext5 ? ext5.toUpperCase() : ext5;
    this.ext3 = ext3 ? ext3.toUpperCase() : ext3;
    this.isDoubleStranded = isDoubleStranded;
    this.isRNA = isRNA;
    this.isCircular = isCircular;
    this.mod_ext3 = mod_ext3 || "";
    this.mod_ext5 = mod_ext5 || "";
  }
}

/**
Creates a new polynucleotide (DNA or RNA) object as JSON
@param {string} sequence - The sequence of the polynucleotide.
@param {string} ext5 - The 5' extension of the current coding strand.
@param {string} ext3 - The 3' extension of the current coding strand.
@param {boolean} isDoubleStranded - Whether the polynucleotide is double stranded.
@param {boolean} isRNA - Whether the polynucleotide is RNA.
@param {boolean} isCircular - Whether the polynucleotide is circular.
@param {string} mod_ext5 - The 5' end modification.
@param {string} mod_ext3 - The 3' end modification.
@return {String} The created polynucleotide object as JSON.
@example
var polynucleotide = polynucleotide("AGCTAGCT", "GATC", "CTAG", true, false, false, null, null);
@customfunction
*/
function polynucleotide(sequence, ext5, ext3, isDoubleStranded, isRNA, isCircular, mod_ext5, mod_ext3) {
  var out = new Polynucleotide(sequence, ext5, ext3, isDoubleStranded, isRNA, isCircular, mod_ext5, mod_ext3);;
  return JSON.stringify(out);
}

/**
* Creates a Polynucleotide object (as JSON) representing a blunt-ended, double-stranded DNA
* such as results from PCR or GBlock synthesis.  It lacks 5' phosphates.
* @param {string} sequence - The coding strand sequence of the dsDNA.
* @return {string} The string representation of the created polynucleotide as JSON.
* @example
* var frag = dsDNA("AGCTAGCT");
* console.log(frag); // Output: '{"sequence":"AGCTAGCT","ext5":null,"ext3":null,"isDoubleStranded":true,
* "isRNA":false,"isCircular":false,"mod_ext5":null,"mod_ext3":null}'
* @customfunction
*/
function dsDNA(sequence) {
  var out = new Polynucleotide(sequence, "", "", true, false, false, "hydroxyl", "hydroxyl");
  return JSON.stringify(out);
}

/**
* Creates a Polynucleotide object (as JSON) representing a linear single stranded DNA (an oligonucleotide)
* @param {string} sequence - The sequence of the polynucleotide oligo.
* @return {string} The string representation of the created polynucleotide oligo.
* @example
* var oligo = oligo("AGCTAGCT");
* console.log(oligo); // Output: '{"sequence":"AGCTAGCT","ext5":null,"ext3":null,"isDoubleStranded":false
* "isRNA":false,"isCircular":false,"mod_ext5":null,"
*/
function oligo(sequence) {
  var out = new Polynucleotide(sequence, null, null, false, false, false, "hydroxyl", null); 
  return JSON.stringify(out);
}

/**
* Creates a Polynucleotide object (as JSON) representing a circular doubkle stranded DNA (a plasmid)
* @param {string} sequence - The sequence of the polynucleotide plasmid.
* @return {string} The string representation of the created polynucleotide plasmid.
* @example
* var plasmid = plasmid("AGCTAGCT");
console.log(plasmid); // Output: '{"sequence":"AGCTAGCT","ext5":null,"ext3":null,"isDoubleStranded":true,"isRNA":false,"isCircular":true,"mod_ext5":null,"mod_ext3":null}'
*/
function plasmid(sequence) {
  var out = new Polynucleotide(sequence, "", "", true, false, true, null, null); 
  return JSON.stringify(out);
}

/**
 * Not to be called from Sheets
 * For resolving a string to a Polynucleotide object
 */
function _resolveToPoly(seqOrJSON) {
  //See if its a JSON already
  try{
      var json = JSON.parse(seqOrJSON);
      return json;
  }
  catch(err) {/*intentionally empty*/}

  //See if its a singular DNA sequence
  if (!_regexDNA.test(seqOrJSON)) {
    throw Error("Cannot resolve " + seqOrJSON);
  }
  
  return new Polynucleotide(seqOrJSON, "", "", true, false, false, "hydroxyl", "hydroxyl");
}

/**
 * Determines whether a given DNA sequence is palindromic.
 * A palindromic DNA sequence is one that reads the same forward and backward when complemented.
 * This function also checks for invalid characters in the input sequence and throws an exception if any are found.
 *
 * @param {string} seq - The input DNA sequence to be checked for palindromicity.
 * @return {boolean} - Returns true if the input DNA sequence is palindromic, false otherwise.
 * @throws {Error} - Throws an error if the input sequence contains characters other than A, T, C, or G.
 *
 * Usage:
 *   const result = isPalindromic("AATT");
 *   console.log(result); // Output: true
 *
 * Example:
 *   1. isPalindromic("AATT") returns true
 *   2. isPalindromic("AGCT") returns false
 *   3. isPalindromic("GAATTC") returns true
 *   4. isPalindromic("AATN") throws an error (invalid character 'N')
 */
function isPalindromic(seq) {
  const complements = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
  };

  // Check for invalid characters and throw an exception if found
  for (const nucleotide of seq) {
    if (!complements.hasOwnProperty(nucleotide)) {
      throw new Error(`Error: Invalid character '${nucleotide}' found in sequence '${seq}'. Sequence must contain only A, T, C, or G.`);
    }
  }

  const reverseComplement = seq.split('').reverse().map(nucleotide => complements[nucleotide]).join('');
  return seq === reverseComplement;
}


/**
 * Calculates the reverse complement of a DNA sequence, including handling of degeneracy codes.
 *
 * @param {string} inseq - The DNA sequence to reverse complement
 * @return {string} The reverse complement of the DNA sequence, or "N/A" if the input contains invalid characters
 */
function revcomp(inseq) {
  if(!inseq.length) {
    return "error on " + inseq;
  }

  var output = "";
  for (let i = inseq.length-1; i >= 0; i--) {
    switch(inseq[i]) {
      case 'A': { output += 'T'; continue; }
      case 'T': { output += 'A'; continue; }
      case 'C': { output += 'G'; continue; }
      case 'G': { output += 'C'; continue; }
      case 'a': { output += 't'; continue; }
      case 't': { output += 'a'; continue; }
      case 'c': { output += 'g'; continue; }
      case 'g': { output += 'c'; continue; }

      case 'B': { output += 'V'; continue; }
      case 'D': { output += 'H'; continue; }
      case 'H': { output += 'D'; continue; }
      case 'K': { output += 'M'; continue; }
      case 'N': { output += 'N'; continue; }
      case 'R': { output += 'Y'; continue; }
      case 'S': { output += 'S'; continue; }
      case 'M': { output += 'K'; continue; }
      case 'V': { output += 'B'; continue; }
      case 'W': { output += 'W'; continue; }
      case 'Y': { output += 'R'; continue; }

      case 'b': { output += 'v'; continue; }
      case 'd': { output += 'h'; continue; }
      case 'h': { output += 'd'; continue; }
      case 'k': { output += 'm'; continue; }
      case 'n': { output += 'n'; continue; }
      case 'r': { output += 'y'; continue; }
      case 's': { output += 's'; continue; }
      case 'm': { output += 'k'; continue; }
      case 'v': { output += 'b'; continue; }
      case 'w': { output += 'w'; continue; }
      case 'y': { output += 'r'; continue; }
      
      default:  throw new Error("Character '" + inseq[i] + "' is not a valid DNA character");
    }
  }
  return output;
}
         

/**
 * Calculates the G/C content of a DNA sequence.
 *
 * @param {string} inseq - The DNA sequence to analyze
 * @return {number} The G/C content of the DNA sequence, between 0 and 1
 */
function gccontent(inseq) {
  inseq = inseq.toUpperCase();
  let gcCount = 0;
  for (let i = 0; i < inseq.length; i++) {
    if (inseq[i] == "G" || inseq[i] == "C") {
      gcCount++;
    }
  }
  return gcCount / inseq.length;
}

/**
 * Calculates the base balance of a DNA sequence.
 *
 * @param {string} inseq - The DNA sequence to analyze
 * @return {number} The base balance of the DNA sequence, between 0 and 1
 */
function basebalance(inseq) {
    inseq = inseq.toUpperCase();
  let baseCounts = {
    A: 0,
    C: 0,
    G: 0,
    T: 0,
  };
  for (let i = 0; i < inseq.length; i++) {
    baseCounts[inseq[i]]++;
  }
  let score = 1;
  for (const base in baseCounts) {
    if (baseCounts[base] == 0) {
      score = 0;
      break;
    }
    score *= baseCounts[base] / inseq.length;
  }
  return 4*Math.pow(score, 1/4);
}

/**
 * Finds the longest streak of repeating bases in a DNA sequence.
 *
 * @param {string} inseq - The DNA sequence to analyze
 * @return {number} The longest streak of repeating bases in the DNA sequence
 */
function maxrepeat(inseq) {
  inseq = inseq.toUpperCase();

  let lastBase = "";
  let streak = 0;
  let maxStreak = 0;
  for (let i = 0; i < inseq.length; i++) {
    if (inseq[i] == lastBase) {
      streak++;
      maxStreak = Math.max(maxStreak, streak);
    } else {
      lastBase = inseq[i];
      streak = 1;
    }
  }
  return maxStreak;
}
