## Check AlphaFold PDB Sequences

Takes target sequences and potential AlphaFold PDBs and records details of sequence matches

### Sequence Checking

Currently limited, carries out following checks:

- Exact sequence match
- Target sequence truncated wrt AlphaFold sequence but exact match of target found _within_ alphafold sequence
  - Residue numbers of overlap are captured
- AlphaFold sequence truncated wrt target sequence but exact match of AlphaFold found _within_ target sequence
  - Residue numbers of overlap are captured

**TODO**

- Fuzzy matching and sub-sequence overlapping for where none of the above matching is found

### Key Functions

- checkAFPDBSequenceForDataFrame
  - Takes a dataframe which must have target sequences and corresponding paths to pdbs
  - Returns sequence match information for each row 
- CheckAlphaFoldPDBSequences_EndToEnd
  - Stitches together dataframe with [uniprotIDs, and sequences] to dataframe with [uniprotIDs, and PDB paths]
    - While this can be used generally it is built with the output from [FetchAlphaFoldPDBs](https://github.com/jackent601/FetchAlphaFoldPDBs) in mind
  - Runs above function on merged dataframe
  - Saves results

### Walkthrough

Explains functions within CheckAFPDBSequences, and tests the demo case.

### Run

Notebook to execute each of the main functions.