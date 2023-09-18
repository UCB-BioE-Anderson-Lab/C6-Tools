/**
 * @file C6-Utils.gs
 * @author J. Christopher Anderson with ChatGPT
 * @copyright 2023 University of California, Berkeley
 * @license See the LICENSE file included in the repository
 * @version 1.0.0
 * @module C6-Utils
 * @description
 * This script provides a collection of utility functions for working with Google Sheets and JSON data.
 * The functions offer various utilities for data manipulation, string concatenation, and JSON parsing.
 */

/**
 * Concatenates an unlimited number of strings, numbers, cells, or arrays using a specified delimiter.
 *
 * @param {...*} args - The values to concatenate. The last argument must be a delimiter string.
 * @returns {string} The concatenated result with the specified delimiter between each value.
 *
 * @example
 * // Assuming the following data in cells A1:A3:
 * // A1: Jenny
 * // A2: Ren
 * // A3: Eva
 * =merge("name",A1:A3,4,", ")
 * // Output: "name, Jenny, Ren, Eva, 4"
 * 
 * @example
 * // Assuming the following DNA sequence in cell B1:
 * // B1: AGCTGACCTGAAGCTAGGTC
 * // Concatenate an EcoRI site (GAATTC) to the 5' end of a 20bp annealing region
 * =merge("GAATTC", B1, "")
 * // Output: "GAATTCAGCTGACCTGAAGCTAGGTC"
 * @customfunction
 */
function merge(...args) {
  if (args.length < 2) {
    throw new Error("At least two arguments are required");
  }

  const delimiter = args[args.length - 1];
  if (typeof delimiter !== "string") {
    throw new Error("The last argument must be a delimiter string");
  }

  const items = args.slice(0, args.length - 1);
  let result = "";

  for (let i = 0; i < items.length; i++) {
    const item = items[i];

    if (Array.isArray(item)) {
      for (let j = 0; j < item.length; j++) {
        const subItem = item[j];
        result += subItem.toString();
        if (j < item.length - 1) {
          result += delimiter;
        }
      }
    } else {
      result += item.toString();
    }

    if (i < items.length - 1) {
      result += delimiter;
    }
  }

  return result;
}

/**
 * Extract a Field from a JSON Object or Array as a String
 * @param {string} objJSON - The JSON object or array expressed as a string
 * @param {string} field - The name of the field to extract from the JSON object
 * @returns {any} The value of the specified field
 * @example
 * field('{"name": "Jenn", "age": 30, "city": "New York"}', "name");
 * returns: "Jenn"
 * @customfunction
*/
function field(objJSON, field) {
  var obj = JSON.parse(objJSON);
  var value = obj[field];
  if (typeof value === "object") {
    return JSON.stringify(value);
  }
  return value;
}

/**
 * Converts a 2D array of string data into a JSON string, representing the contents of the array data in JSON format.
 * The function assumes that the input array follows a specific structure, with field names ending with a colon (:).
 * It processes objects, arrays, and single values, handling 'null' and 'undefined' values accordingly.
 *
 * Usage example:
 *
 * var inputArray = [
 *   ["oligos:", "", "name:", "sequence:", "description:", "supplier:", "date:"],
 *   ["", "G00101", "ATTACCGCCTTTGAGTGAGC", "BBa_G00101 sequencing primer", "IDT", "12/23/23"],
 *   ["", "cyz4", "CACGATCAGTCCGCGTTTG", "Construction of pCyz4", "", ""],
 *   ["supplier:", "IDT", "", "", "", "", ""],
 *   ["date:", "12/23/23", "", "", "", "", ""]
 * ];
 * var jsonString = makeJSON(inputArray);
 * Logger.log(jsonString);
 *
 * Output:
 *
 * {"oligos":[{"name":"G00101","sequence":"ATTACCGCCTTTGAGTGAGC","description":"BBa_G00101 sequencing primer"},{"name":"cyz4","sequence":"CACGATCAGTCCGCGTTTG","description":"Construction of pCyz4"}],"supplier":"IDT","date":"12/23/23"}
 *
 * @param {Array<Array<string>>} inputArray The 2D array of strings representing the data to be converted to JSON.
 * @return {string} A JSON string representing the contents of the input data.
 * @customfunction
 */
function makeJSON(inputArray) {
  // Convert inputArray to a new array with .toString() values
  var stringArray = inputArray.map(function (row) {
    return row.map(function (value) {
      return value.toString();
    });
  });

  function analyzeStructure(rangeArray) {
    var structure = {};

    for (var i = 0; i < rangeArray.length; i++) {
      var key = rangeArray[i][0];
      if (key.endsWith(':')) {
        key = key.slice(0, -1);
        if (rangeArray[i][1] === '') {
          structure[key] = 'array';
          i += 1; // Skip header row
        } else {
          structure[key] = 'value';
        }
      }
    }

    return structure;
  }

  function processRange(rangeArray, structure) {
    var resultObj = {};
    var i = 0;

    while (i < rangeArray.length) {
      var row = rangeArray[i];
      var key = row[0];

      if (key.endsWith(':')) {
        key = key.slice(0, -1);

        if (structure[key] === 'value') {
          resultObj[key] = row[1] === 'null' ? null : row[1] === 'undefined' ? undefined : row[1];
          i++;
        } else if (structure[key] === 'array') {
          var arrayKey = key;
          var arrayData = [];
          var headerRow = rangeArray[i + 1];
          i += 2;

          while (i < rangeArray.length && !rangeArray[i][0].endsWith(':')) {
            var obj = {};
            for (var k = 1; k < headerRow.length; k++) {
              var fieldKey = headerRow[k].endsWith(':') ? headerRow[k].slice(0, -1) : headerRow[k];
              if (rangeArray[i][k] !== undefined) {
                obj[fieldKey] = rangeArray[i][k] === 'null' ? null : rangeArray[i][k] === 'undefined' ? undefined : rangeArray[i][k];
              }
            }
            arrayData.push(obj);
            i++;
          }

          resultObj[arrayKey] = arrayData;
        }
      } else {
        i++;
      }
    }

    return resultObj;
  }

  var structure = analyzeStructure(stringArray);
  return JSON.stringify(processRange(stringArray, structure));
}

/**
 * Parses a JSON string into a 2D array of strings, representing the contents of the JSON data in tabular form.
 * Nested arrays and objects are flattened into columns and indented to show hierarchy.
 * Null and undefined values are replaced with the strings 'null' or 'undefined'.
 *
 * Usage example:
 *
 * var inputString = '{"oligos":[{"name":"G00101","sequence":"ATTACCGCCTTTGAGTGAGC","description":"BBa_G00101 sequencing primer"},{"name":"cyz4","sequence":"CACGATCAGTCCGCGTTTG","description":"Construction of pCyz4"}], "supplier":"IDT", "date":"12/23/23"}';
 * var outputArray = parseJSON(inputString);
 * Logger.log(outputArray);
 *
 * Output:
 *
 * [
 *   ["oligos:", "", "", ""],
 *   ["", "name:", "sequence:", "description:"],
 *   ["", "G00101", "ATTACCGCCTTTGAGTGAGC", "BBa_G00101 sequencing primer"],
 *   ["", "cyz4", "CACGATCAGTCCGCGTTTG", "Construction of pCyz4"],
 *   ["supplier:", "IDT", "", ""],
 *   ["date:", "12/23/23", "", ""]
 * ]
 *
 * @param {string} inputString The JSON string to parse.
 * @return {Array<Array<string>>} A 2D array of strings representing the contents of the JSON data in tabular form.
 * @customfunction
 */
function parseJSON(inputString) {
  var jsonData = JSON.parse(inputString);
  var outputArray = [];

  function flattenObject(obj, prefix) {
    prefix = prefix || '';
    Object.keys(obj).forEach(function(key) {
      var fullKey = prefix + key;
      if (typeof obj[key] === 'object' && !Array.isArray(obj[key]) && obj[key] !== null) {
        outputArray.push([fullKey + ':']);
        flattenObject(obj[key], fullKey + '.');
      } else if (Array.isArray(obj[key])) {
        outputArray.push([fullKey + ':']);
        var keys = Object.keys(obj[key][0] || {});
        outputArray.push(['', ...keys.map(function(k) { return k + ':'; })]);
        for (var i = 0; i < obj[key].length; i++) {
          var values = Object.values(obj[key][i] || {}).map(function(v) { return v === null || v === undefined ? v + '' : v; });
          outputArray.push(['', ...values]);
        }
      } else {
        var value = obj[key] === null || obj[key] === undefined ? obj[key] + '' : obj[key];
        outputArray.push([fullKey + ':', value]);
      }
    });
  }
  
  if (typeof jsonData === 'object' && jsonData !== null) {
    flattenObject(jsonData);
  } else {
    throw new Error('Invalid input: Not a JSON object.');
  }
  
  return outputArray;
}


