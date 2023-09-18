/**
 * @file C6-Sim.gs
 * @author J. Christopher Anderson with ChatGPT
 * @copyright 2023 University of California, Berkeley
 * @license See the LICENSE file included in the repository
 * @version 1.0.0
 * @module C6-Sim
 * @description
 * This script provides a collection of functions for expressing construction files and simulating
 * molecular biology operations including PCR and assembly reactions.
 *
 * @requires C6-Utils
 * @requires C6-Seq
 * @requires C6-Oligos
 */

/**
 * A Construction File (CF) is a structured format for specifying a series of molecular biology construction steps, 
 * such as PCR, assembly, digestion, ligation, and transformation. It is designed to facilitate communication between 
 * researchers and computer programs, allowing users to simulate or perform complex DNA manipulations.
 *
 * In this project, a CF is expressed as JSON. There is a function parseCF which can read in data and convert it to
 * this JSON format. There is another function, simCF, which can input this JSON and simulate the steps. You can
 * also simulate individual steps one at a time by invoking the PCR, assemble, etc. functions.
 *
 * Usage:
 * CFs are used to plan, simulate, and document molecular biology experiments. Researchers can design and share
 * their construction steps in a standardized format, and software programs can parse, simulate, and visualize the
 * planned steps, providing an efficient way to manage and analyze experimental data.
 *
 * Syntax:
 * A Construction File is represented as a JSON object, containing two main elements: 'steps'
 * and 'sequences'. The 'steps' is an array of objects, where each object represents a construction
 * step with its associated operation, input sequences, and output product. The 'sequences' is an object
 * containing key-value pairs, where each key is a unique identifier for a DNA sequence, and the value is the
 * actual sequence.
 *
 * Example:
 * Here's a simple example of a Construction File that demonstrates PCR and assembly steps.
 *
 * {
 *   "steps": [
 *     {
 *       "operation": "PCR",
 *       "output": "P6",
 *       "forward_oligo": "P6libF",
 *       "reverse_oligo": "P6libR",
 *       "template": "pTP1",
 *       "product_size": 3583
 *     },
 *     {
 *       "operation": "Assemble",
 *       "output": "pP6",
 *       "dnas": ["P6"],
 *       "enzyme": "BsaI"
 *     }
 *   ],
 *   "sequences": {
 *     "P6libF": "ccaaaggtctcATTATANNNNNNNNNNNNNNNNNTGTCAANNNNGAacccaggactcctcgaagtcgttcttaagacaac",
 *     "P6libR": "cagttGGTCTCAATAATNNNNNNANNNNGTtagtatttctcctcgtctacggttaactgatactc",
 *     "pTP1": "ATTACCGCCTTTGAGTGG"
 *   }
 * }
 *
 * In this example, the 'steps' array has two steps: PCR and Assemble. The PCR step uses forward and
 * reverse oligos "P6libF" and "P6libR", with "pTP1" as the template. The PCR product is named "P6". The Assemble
 * step uses the "P6" PCR product and the "BsaI" enzyme to create a final output named "pP6". The 'sequences'
 * object contains the sequences for "P6libF", "P6libR", and "pTP1".
 *
 * @typedef {Object} ConstructionFile
 * @property {Array.<PCR|Assemble|GoldenGate|Gibson|Transform|Digest|Ligate>} steps - An array of construction steps, where each step is an operation object.
 * @property {Object} sequences - An object containing key-value pairs of sequence names and their corresponding DNA sequences.
 *
 * @typedef {Object} PCR
 * @property {'PCR'} operation - The type of operation.
 * @property {string} output - The output product of the operation.
 * @property {string} forward_oligo - Forward primer used in PCR operation.
 * @property {string} reverse_oligo - Reverse primer used in PCR operation.
 * @property {string} template - The template DNA used in PCR operation.
 * @property {number} product_size - The expected product size in PCR operation.
 *
 * @typedef {Object} Assemble
 * @property {'Assemble'} operation - The type of operation.
 * @property {string} output - The output product of the operation.
 * @property {Array.<string>} dnas - An array of DNA parts used in the Assemble operation.
 * @property {string} enzyme - The enzyme used in the Assemble operation.
 *
 * @typedef {Object} Transform
 * @property {'Transform'} operation - The type of operation.
 * @property {string} output - The output product of the operation.
 * @property {string} dna - The DNA used in the Transform operation.
 * @property {string} strain - The bacterial strain used in the Transform operation.
 * @property {string} antibiotics - The antibiotics used in the Transform operation.
 * @property {number} [temperature] - The temperature used in the Transform operation (optional).
 *
 * @typedef {Object} Digest
 * @property {'Digest'} operation - The type of operation.
 * @property {string} output - The output product of the operation.
 * @property {string} dna - The DNA used in the Digest operation.
 * @property {number} fragSelect - The index, counted from zero, of the output fragment
 * @property {Array.<string>} enzymes - The enzymes used in the Digest operation.
 *
 * @typedef {Object} Ligate
 * @property {'Ligate'} operation - The type of operation.
 * @property {string} output - The output product of the operation.
 * @property {Array.<string>} dnas - An array of DNA parts used in the Ligate operation.
 */

/**
 * parseCF - A function to parse construction and sequence data from various input formats.
 * 
 * Usage:
 * 
 * const output = parseCF(...blobs);
 * 
 * Arguments:
 * 
 * - blobs: One or more inputs containing construction and sequence data. Each input can be:
 *   - A single cell value (string)
 *   - A 1D array of strings (e.g., a row or column of cell values)
 *   - A 2D array of strings (e.g., a range of cell values)
 *   
 * The function processes the input data, identifies construction and sequence data,
 * and outputs a JSON string containing the parsed data organized into steps
 * and sequences.
 * 
 * Example input data formats:
 * 
 * 1. Single cell value:
 * 
 * "PCR P6libF P6libR on pTP1, P6"
 * 
 * 2. 1D array (e.g., row of cell values):
 * 
 * ["PCR", "P6libF", "P6libR", "on", "pTP1", "P6"]
 * 
 * 3. 2D array (e.g., range of cell values):
 * 
 * [
 *   ["PCR", "P6libF", "P6libR", "on", "pTP1", "P6"],
 *   ["Assemble", "pTP1", "P6", "pP6", "P6"]
 * ]
 * 
 * Example output:
 * 
 * {
 *   "steps": [
 *     {
 *       "operation": "PCR",
 *       "output": "P6",
 *       "forward_oligo": "P6libF",
 *       "reverse_oligo": "P6libR",
 *       "template": "pTP1"
 *     },
 *     {
 *       "operation": "Assemble",
 *       "output": "pP6",
 *       "dnas": ["pTP1", "P6"],
 *       "enzyme": "P6"
 *     }
 *   ],
 *   "sequences": {}
 * }
 * 
 */
function parseCF(...blobs) {
    const operationRegex = /^(PCR|pcr|digest|Digest|ligate|Ligate|assemble|Assemble|gibson|Gibson|goldengate|GoldenGate|blunt|Blunt|transform|Transform)/;
    const sequenceDataRegex = /^[ACGTRYSWKMBDHVNUacgtryswkmbdhvnu*]+$/;

    // Helper function to tokenize construction data
    function preprocessData(data) {
        if (Array.isArray(data)) {
            if (data.every(item => Array.isArray(item))) {
                // Handle 2D array
                return data.map(row => row.map(cell => cell.toString()).join('\t')).join('\n');
            } else if (data.every(item => typeof item === "string" || typeof item === "number")) {
                // Handle 1D array
                return data.map(cell => cell.toString()).join('\t');
            } else {
                throw new Error("Unsupported input type for preprocessData function");
            }
        } else {
            // Handle single-cell value
            return data.toString();
        }
    }

    function tokenize(text) {
        // Split the text into tokens using whitespace, commas, and slashes as separators
        let tokens = text.split(/[\s,/()]+/);
        // Remove any punctuation from each token
        tokens = tokens.map(token => token.replace(/[(),]/g, ""));
        // Remove the words "on" and "with" from the tokens
        tokens = tokens.filter(token => !["on", "with"].includes(token));
        // Remove any empty tokens resulting from punctuation removal
        tokens = tokens.filter(token => token !== "");
        return tokens;
    }

    // Preprocess input data into a single text string
    let singleblob = "";
    for (const blob of blobs) {
        singleblob += preprocessData(blob) + '\n';
    }

    // Tokenize the singleblob into a 2D array of strings
    const preprocessedData = singleblob.trim().split('\n').map(line => tokenize(line));

    // Helper function to parse construction data
    function parseConstructionData(data) {
        function processTokens(tokens) {
            const operation = tokens[0];

            let step = {
                operation
            };

            console.log("Parsing step: " + operation);

            switch (operation.toLowerCase()) {
                case "pcr":
                    step.operation = 'PCR';
                    step.output = tokens[tokens.length - 1];
                    step.forward_oligo = tokens[1];
                    step.reverse_oligo = tokens[2];
                    step.template = tokens[3];
                    if (tokens[4]) step.product_size = parseInt(tokens[4], 10);
                    break;
                case "assemble":
                    step.operation = 'Assemble';
                    step.output = tokens[tokens.length - 1];
                    step.dnas = tokens.slice(1, tokens.length - 2);
                    step.enzyme = tokens[tokens.length - 2];
                    break;
                case "goldengate":
                    step.operation = 'Assemble';
                    step.output = tokens[tokens.length - 1];
                    step.dnas = tokens.slice(1, tokens.length - 2);
                    step.enzyme = tokens[tokens.length - 2];
                    console.log("GoldenGate triggered");
                    break;
                case "gibson":
                    step.operation = 'Assemble';
                    step.output = tokens[tokens.length - 1];
                    step.dnas = tokens.slice(1, tokens.length - 1);
                    step.enzyme = "Gibson";
                    break;
                case "ligate":
                    step.operation = 'Ligate';
                    step.output = tokens[tokens.length - 1];
                    step.dnas = tokens.slice(1, tokens.length - 1);
                    break;
                case "digest":
                    step.operation = 'Digest';
                    step.output = tokens[tokens.length - 1];
                    step.dna = tokens[1];
                    step.enzymes = tokens[2];
                    step.fragselect = parseInt(tokens[3], 10);
                    break;
                case "transform":
                    step.operation = 'Transform';
                    step.output = tokens[4];
                    step.dna = tokens[1];
                    step.strain = tokens[2];
                    step.antibiotics = tokens[3];
                    step.temperature = parseFloat(tokens[5]);
                    break;
            }

            return step;
        }

        // Check if the input is a range of cells
        if (data instanceof Array && data.length > 0 && data[0].length > 1) {
            return data.map(row => processTokens(row));
        }
        // Check if the input is a series of lines
        if (data instanceof Array && data.length > 1 && data[0].length === 1) {
            // Split each element of the row into tokens and filter out empty tokens
            return data.map(line => processTokens(tokenize(line).filter(token => token !== "")));
        }
        // Check if the input is a string
        if (typeof data === "string") {
            // Split the string into lines, remove empty lines and comments
            const lines = data.split(/\r?\n/).filter(line => !line.match(/^\s*$|\/\//));
            // Split each line into an array of tokens and filter out empty tokens
            return lines.map(line => processTokens(tokenize(line).filter(token => token !== "")));
        }
        // If none of the above cases match, throw an error
        throw new Error("Unable to parse construction data");
    }

    // Helper function to validate the input data
    function validateData(steps, sequences) {
        // Add your implementation here
    }

    // Initialize construction and sequence data storage
    let steps = [];
    let sequences = {};

    // Iterate over preprocessedData
    for (const row of preprocessedData) {
        if (operationRegex.test(row[0])) {
            const parsedConstructionData = parseConstructionData([row]);
            steps.push(...parsedConstructionData);
        } else {
            const name = row[0];
            const sequence = row.slice(1).join('');
            if (sequenceDataRegex.test(sequence)) {
                sequences[name] = sequence;
            }
        }
    }

    // Validate the input data
    validateData(steps, sequences);

    // Bundle all the info together and output an organized object or array
    const output = {
        steps,
        sequences
    };

    return JSON.stringify(output);
}

/**
 * PCR function predicts the sequence of a PCR product by inputting forward oligo sequence, reverse oligo sequence, and template sequence.
 *
 * The algorithm assumes that the last 18 bp on the 3' end of the oligos exactly match the template, but the 5' end of the oligos may not match.
 * The function first identifies where the forward oligo will anneal to the template by looking for the 18 bp match.
 * It then rotates the template sequence such that it begins with the annealing region of the forward oligo.
 * Then it looks for the annealing site of the reverse oligo by invoking the function revcomp(sequence) which inputs the reverse oligo sequence 
 * and outputs the reverse complement. 
 * The first 18 bp of that revcomp sequence should match the template where it will anneal.
 * Based on the indices of the annealing sites, the final PCR product is calculated from the entire forward sequence, the region between the 
 * annealing regions on the rotated template, and the entire reverse complement of the reverse oligo
 *
 * @param {string} forwardSeq - The forward oligo sequence.
 * @param {string} reverseSeq - The reverse oligo sequence.
 * @param {string} templateSeq - The template sequence.
 *
 * @returns {string} finalProduct - The predicted PCR product.
 */
function PCR(forwardSeq, reverseSeq, templateSeq) {
    try {
        forwardSeq = resolveToSeq(forwardSeq);
    } catch (err) {
        throw new Error('PCR unable to parse forward primer');
    }
    try {
        reverseSeq = resolveToSeq(reverseSeq);
    } catch (err) {
        throw new Error('PCR unable to parse reverse primer');
    }
    try {
        templateSeq = resolveToSeq(templateSeq);
    } catch (err) {
        throw new Error('PCR unable to parse template sequence');
    }

    // Find index of 18 bp match on 3' end of forward oligo and template
    var foranneal = forwardSeq.slice(-18);
    // console.log(foranneal);

    var forwardMatchIndex = templateSeq.indexOf(foranneal);
    // console.log(forwardMatchIndex);
    if (forwardMatchIndex === -1) {
        templateSeq = revcomp(templateSeq);
        forwardMatchIndex = templateSeq.indexOf(foranneal);
        if (forwardMatchIndex === -1) {
            throw new Error("Forward oligo does not exactly anneal to the template")
        }
    }

    // Rotate template sequence to begin with annealing region of forward oligo
    var rotatedTemplate = templateSeq.slice(forwardMatchIndex) + templateSeq.slice(0, forwardMatchIndex);

    // Find reverse complement of reverse oligo
    var reverseComp = revcomp(reverseSeq);

    // Find index of 18 bp match on 3' end of reverse complement and rotated template
    var revanneal = reverseComp.slice(0, 18);
    var reverseMatchIndex = rotatedTemplate.indexOf(revanneal);
    if (reverseMatchIndex === -1) {
        throw new Error("Reverse oligo does not exactly anneal to the template")
    }

    // Concatenate entire forward oligo, region between annealing regions on rotated template, and entire 
    // reverse complement of reverse oligo to obtain final PCR product
    var finalProduct = forwardSeq + rotatedTemplate.slice(18, reverseMatchIndex) + reverseComp;

    return finalProduct;
}

/**
 * The following blocks describe the cutting pattern of commonly
 * used restriction enzymes.  They are used in the Assemble and
 * Digest simulations, then also during silent site removal (removeSites)
 */
const restrictionEnzymes = {
    AarI: {
        recognitionSequence: "CACCTGC",
        cut5: 4,
        cut3: 8
    },
    BbsI: {
        recognitionSequence: "GAAGAC",
        cut5: 2,
        cut3: 6
    },
    BsaI: {
        recognitionSequence: "GGTCTC",
        cut5: 1,
        cut3: 5
    },
    BsmBI: {
        recognitionSequence: "CGTCTC",
        cut5: 1,
        cut3: 5
    },
    SapI: {
        recognitionSequence: "GCTCTTC",
        cut5: 1,
        cut3: 4
    },
    BseRI: {
        recognitionSequence: "GAGGAG",
        cut5: 10,
        cut3: 8
    },
    BamHI: {
        recognitionSequence: "GGATCC",
        cut5: -5,
        cut3: -1
    },
    BglII: {
        recognitionSequence: "AGATCT",
        cut5: -5,
        cut3: -1
    },
    EcoRI: {
        recognitionSequence: "GAATTC",
        cut5: -5,
        cut3: -1
    },
    XhoI: {
        recognitionSequence: "CTCGAG",
        cut5: -5,
        cut3: -1
    },
    SpeI: {
        recognitionSequence: "ACTAGT",
        cut5: -5,
        cut3: -1
    },
    XbaI: {
        recognitionSequence: "TCTAGA",
        cut5: -5,
        cut3: -1
    },
    PstI: {
        recognitionSequence: "CTGCAG",
        cut5: -1,
        cut3: -5
    },
}

//Calculate the reverse complement of the recognition sequence as well as
//whether the sticky end generated is a 5' extension (true) or 3' (false)
for (enzName in restrictionEnzymes) {
    const enzyme = restrictionEnzymes[enzName];
    enzyme.recognitionRC = revcomp(enzyme.recognitionSequence);
    enzyme.isFivePrime = enzyme.cut5 < enzyme.cut3;
}

/**
 * Assembles a set of DNA sequences using the Gibson or Golden Gate Assembly method
 *
 * @param {Array<string|Array<string>>} dnaBlobs An array of DNA sequences or sequence fragments to be assembled
 * @param {string} enzyme The restriction enzyme used to cut the DNA sequences, or 'gibson' for Gibson assembly
 * @return {string} The final assembled DNA sequence
 *
 * Approach:
 * 1. Verify that the input enzyme is a valid restriction enzyme
 * 2. Use the restriction enzyme to cut the DNA sequences and collect the digestion fragments
 * 3. Sort the digestion fragments based on their sticky end compatibility
 * 4. Check that each sticky end has a matching sticky end in the next fragment
 * 5. Concatenate the cut fragments to form the final assembled DNA sequence
 */
function assemble(...dnaBlobs) {
    // dnaBlobs = ["CCATAGGTCTCAGCTTCTACTAGAGCATAAGCGTGGCTTAACAATTCCCTACTAGAGACCTTGTC","CCATAGGTCTCATACTATTAAGGTGGAGAAAGGTCAGGCCGGCTTAGAGACCTTGTC","BsaI"];
    const enzyme = dnaBlobs.pop();

    // Flatten the array of dnaBlobs into a single array of cleaned up DNA sequences
    const dnaSequences = dnaBlobs.flat().map(sequence => {
        if (Array.isArray(sequence)) {
            // Clean up each sequence in the array
            return sequence.map(subsequence => resolveToSeq(subsequence));
        } else {
            // Clean up the single sequence
            return resolveToSeq(sequence);
        }
    }).flat();

    // check if the enzyme is in the restrictionEnzymes object, or relay to Gibson
    if (restrictionEnzymes.hasOwnProperty(enzyme)) {
        return goldengate(dnaSequences, enzyme);
    } else {
        return gibson(dnaSequences);
    }
}

/**
 * Assembles DNA sequences using the Golden Gate assembly method.
 *
 * @param {Array} seqs - an array of DNA sequences to be assembled, supplied as strings.
 * @param {string} enzyme The restriction enzyme used to cut the DNA sequences
 *
 * @returns {string} - the assembled DNA sequence.
 *
 * @throws {Error} - if the input is not a non-empty array of DNA sequences, or if the assembly does not resolve to
 *                   a single product, or if the products do not assemble correctly.
 */
function goldengate(dnaSequences, enzyme) {

    // check if the enzyme is in the restrictionEnzymes object, or relay to Gibson
    if (!restrictionEnzymes.hasOwnProperty(enzyme)) {
        throw new Error(enzyme + " is not a supported Golden Gate enzyme");
    }

    // Initialize variables to store the restriction enzyme sites and the assembly product
    var assemblyProduct = "";
    var digestionFragments = [];

    // Get the restriction enzyme details
    var enzymeDetails = restrictionEnzymes[enzyme];
    var restrictionSequence = enzymeDetails.recognitionSequence;
    var revRestrictionSequence = enzymeDetails.recognitionRC;
    var cut5 = enzymeDetails.cut5;
    var cut3 = enzymeDetails.cut3;

    // Loop through the DNA sequences and collect digestion fragments
    dnaSequences.forEach(sequence => {
        sequence = resolveToSeq(sequence);

        // Find the restriction enzyme sites in the current sequence
        const enzymeSites = sequence.split(restrictionSequence).length - 1;
        const revEnzymeSites = sequence.split(revRestrictionSequence).length - 1;

        // Declare enzymeSite and revEnzymeSite
        const enzymeSite = sequence.indexOf(restrictionSequence);
        const revEnzymeSite = sequence.indexOf(revRestrictionSequence);

        // If the enzyme is not found in the sequence, throw an error
        if (enzymeSites === 0) {
            throw new Error(`Enzyme site ${restrictionSequence} not found in sequence ${sequence}`);
        }
        if (revEnzymeSites === 0) {
            throw new Error(`Reverse Enzyme site ${revRestrictionSequence} not found in sequence ${sequence}`);
        }

        // If there is more than one forward or reverse site in the sequence, throw an error
        if (enzymeSites > 1) {
            throw new Error(`More than one forward enzyme site ${restrictionSequence} found in sequence ${sequence}`);
        }
        if (revEnzymeSites > 1) {
            throw new Error(`More than one reverse enzyme site ${revRestrictionSequence} found in sequence ${sequence}`);
        }

        // Confirm that the forward site comes before the reverse complement site
        if (revEnzymeSite < enzymeSite) {
            throw new Error(`Reverse enzyme site found before forward enzyme site in sequence ${sequence}`);
        }

        // Extract the cut fragment, stickyEnd5 and stickyEnd3
        let cutFragment, stickyEnd5, stickyEnd3;
        if (enzymeDetails.isFivePrime) {
            cutFragment = sequence.substring(enzymeSite + restrictionSequence.length + cut3, revEnzymeSite - cut3);
            stickyEnd5 = sequence.substring(enzymeSite + restrictionSequence.length + cut5, enzymeSite + restrictionSequence.length + cut3);
            stickyEnd3 = sequence.substring(revEnzymeSite - cut3, revEnzymeSite - cut5);
        } else {
            cutFragment = sequence.substring(enzymeSite + restrictionSequence.length + cut5, revEnzymeSite - cut5);
            stickyEnd5 = sequence.substring(enzymeSite + restrictionSequence.length + cut3, enzymeSite + restrictionSequence.length + cut5);
            stickyEnd3 = sequence.substring(revEnzymeSite - cut5, revEnzymeSite - cut3);
        }

        digestionFragments.push({
            fragment: cutFragment,
            stickyEnd5: stickyEnd5,
            stickyEnd3: stickyEnd3
        });
    });


    // Sort the digestion fragments based on the sticky ends
    digestionFragments.sort((a, b) => {
        if (a.stickyEnd5 === b.stickyEnd3) {
            return 0;
        } else if (a.stickyEnd5 < b.stickyEnd3) {
            return -1;
        } else {
            return 1;
        }
    });

    //Validate that all sticky ends are non-palindromic
    digestionFragments.forEach(fragment => {
        if (isPalindromic(fragment.stickyEnd5) || isPalindromic(fragment.stickyEnd3)) {
            throw new Error(`Palindromic sticky ends found in fragment ${fragment.fragment}`);
        }
    });

    // Check if there is a way to assemble the fragments without including all the fragments
    if (digestionFragments.length > 1) {
        const stickyEndCounts = {};

        digestionFragments.forEach(fragment => {
            if (!stickyEndCounts.hasOwnProperty(fragment.stickyEnd5)) {
                stickyEndCounts[fragment.stickyEnd5] = {
                    count5: 0,
                    count3: 0
                };
            }
            if (!stickyEndCounts.hasOwnProperty(fragment.stickyEnd3)) {
                stickyEndCounts[fragment.stickyEnd3] = {
                    count5: 0,
                    count3: 0
                };
            }
            stickyEndCounts[fragment.stickyEnd5].count5++;
            stickyEndCounts[fragment.stickyEnd3].count3++;
        });

        for (const stickyEnd in stickyEndCounts) {
            if (stickyEndCounts[stickyEnd].count5 > 1 || stickyEndCounts[stickyEnd].count3 > 1) {
                throw new Error('Some fragments have the same sticky ends, which can lead to incorrect assemblies');
            }
        }
    }

    // Validate that the sticky ends match between fragments
    for (var i = 0; i < digestionFragments.length - 1; i++) {
        if (digestionFragments[i].stickyEnd3 !== digestionFragments[i + 1].stickyEnd5) {
            throw new Error(`Sticky ends do not match between fragments 
              ${digestionFragments[i].fragment} and ${digestionFragments[i + 1].fragment}`);
        }
    }
    if (digestionFragments[0].stickyEnd5 !== digestionFragments[digestionFragments.length - 1].stickyEnd3) {
        throw new Error(`Sticky ends do not match between first and last fragments 
          ${digestionFragments[0].fragment} and ${digestionFragments[digestionFragments.length - 1].fragment}`);
    }

    //Assemble the final sequence and return
    var finalSeq = "";
    for (var i = 0; i < digestionFragments.length; i++) {
        finalSeq += digestionFragments[i].stickyEnd5;
        finalSeq += digestionFragments[i].fragment;
    }
    return finalSeq;
}

/**
 * Assembles DNA sequences using the Gibson assembly method.
 * It is also the default algorithm for 'assemble' function.
 * It is also appropriate for SOEing and yeast assembly predictions.
 *
 * @param {Array} seqs - an array of DNA sequences to be assembled, supplied as strings.
 * @param {boolean} check_circular - whether to check if the assembled product is circular. If set to `false`,
 *                                   the function will not check if the product is circular and will return a linear
 *                                   product. Defaults to `true`.
 *
 * @returns {string} - the assembled DNA sequence.
 *
 * @throws {Error} - if the input is not a non-empty array of DNA sequences, or if the assembly does not resolve to
 *                   a single product, or if the products do not assemble correctly, or if the assembled product is
 *                   not circular and `check_circular` is set to `true`.
 */
const HOMOLOGY_LENGTH = 20; //20 bp homology arms required
function gibson(seqs, check_circular) {
    if (!Array.isArray(seqs) || seqs.length === 0) {
        throw new Error("Invalid input: expected non-empty array of DNA sequences");
    }

    //Do sequence cleanup and checks
    for (var i in seqs) {
        seqs[i] = resolveToSeq(seqs[i]);
    }

    //Default to requiring a circular product
    var checkcirc = true;
    if (check_circular === false) {
        checkcirc = false;
    }

    //Do iterations of assembly, one pair of fragments at a time
    var startList = [...seqs];
    var endList = [];
    var isCircular = false;

    //Iterate until it converges on a single sequence
    while (startList.length > 1) {
        //Compare the sequences for homology pairwise
        for (seq1 of startList) {
            const threePrime = seq1.substring(seq1.length - HOMOLOGY_LENGTH);
            for (seq2 of startList) {
                if (seq2 === seq1) {
                    continue;
                }
                var startIndex = seq2.indexOf(threePrime);
                if (startIndex === -1) {
                    continue;
                }

                //It found a match; add the extended product to the list for the next cycle
                var newseq = seq1 + seq2.substring(startIndex + threePrime.length);
                endList.push(newseq);
            }
        }

        //Reset the startList to be the sequences extended in previous cycle
        //If they are the same length, it is a circular product
        if (endList.length === startList.length) {
            startList = endList.slice(0, -1);
            endList = [];
            isCircular = true;

            //If it were linear, or a later cycle, this would be the case
        } else if (endList.length < startList.length) {
            startList = endList;
            endList = [];

            //If there are more seqs in endList than startList, then it is not converging
        } else {
            throw new Error("Products do not assembly correctly, multiply assembly junctions present");
        }
    }

    if (startList.length != 1) {
        throw new Error("Gibson assembly did not resolve to a single product");
    }

    var outseq = startList[0];

    //If it was detected as circular earlier, recircularize with one last recombination
    if (isCircular) {
        const threePrime = outseq.substring(outseq.length - HOMOLOGY_LENGTH);
        var startIndex = outseq.indexOf(threePrime);
        outseq = outseq.substring(startIndex + HOMOLOGY_LENGTH);

        //Otherwise, by default throw an error unless explicitly requested not to
    } else {
        if (checkcirc) {
            throw new Error("Products do not assemble into a circular product.  If you are expecting a linear product, pass in check_circular = FALSE");
        }
    }

    return outseq;
}

/**
simCF - A function that simulates a series of molecular biology construction steps given a JSON string of a construction file.

Each row represents a product with the first column being the name and the second column being the sequence.
The function first parses the input JSON string and retrieves the construction steps and sequence data. It then iterates through
the construction steps in order, simulating each step with the appropriate simulator function. The product sequences for each step
are collected in an array. Finally, the collected product sequences are converted into a 2D string array with the left column being
the product name and the right column being the sequence.
Note: This function relies on the existence of simulator functions (PCR, assemble, etc.) in the same scope.

@param {string} jsonString - A JSON string of a construction file in the format outputted from parseCF.
@returns {Array<Array<string>>} outputTable - A 2D string array containing the product names and their sequences.
*/
function simCF(jsonString) {
    const inputData = JSON.parse(jsonString);
    const steps = inputData.steps;
    const sequences = inputData.sequences;
    const products = [];

    if (!sequences || Object.keys(sequences).length === 0) {
        return "Error: Sequence data is missing. Please include sequence data in the input JSON.";
    }

    function lookupSequence(key) {
        const foundProduct = products.find((product) => product.name === key);
        if (foundProduct) {
            return foundProduct.sequence;
        }

        const foundSequence = sequences[key];
        if (foundSequence) {
            return foundSequence;
        }

        throw new Error(`Missing sequence for key: ${key}`);
    }

    for (let i = 0; i < steps.length; i++) {
        const step = steps[i];

        switch (step.operation) {
            case 'PCR': {
                const forwardOligoSeq = lookupSequence(step.forward_oligo);
                const reverseOligoSeq = lookupSequence(step.reverse_oligo);
                const templateSeq = lookupSequence(step.template);

                const product = PCR(forwardOligoSeq, reverseOligoSeq, templateSeq);
                products.push({
                    name: step.output,
                    sequence: product
                });
            }
            break;

            case 'Assemble': {
                const dnaSequences = step.dnas.map((dnaKey) => lookupSequence(dnaKey));

                const product = assemble(dnaSequences, step.enzyme);
                products.push({
                    name: step.output,
                    sequence: product
                });
            }
            break;

            case 'Transform': {
                const dnaSequence = lookupSequence(step.dna);
                const product = dnaSequence;
                products.push({
                    name: step.output,
                    sequence: product
                });
            }
            break;

            // ... add more cases for other operations as needed

            default:
                throw new Error(`Unsupported operation: ${step.operation}`);
        }
    }

    const outputTable = products.map((product) => [product.name, product.sequence]);
    return outputTable;
}