import pandas as pd
import os
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1 # convert 3 letter AA seq to 1 letter AA seq

# ================================================================================================================================
#   Checking PDB Sequences
# ================================================================================================================================

def get1LetterPrimarySequenceFromModel(model):
    """Requires seq1 from bio.SeqUtils, model from bio.PDB.PDBParser"""
    return seq1("".join([r.get_resname() for r in list(model.get_residues())]))

def checkAFPDBSequence(sourceSequence, AF_PDB_obj, debug = True):
    """    
    Checks if source sequence matches or a truncated matches AlphaFold PDB sequence. In case of truncated match finds the sequence numbers of overlap
    Currently only detects exact matches, need to implement fuzzy matching incase there are still overlapping matches between source and AF

    Could probably be programmed in fewer lines with clever flags to return dictionary at once but this is straightforward 
    Returns a dictionary with information on the match found
    """
    # Get primary sequence to check for match (AlphaFold only has a single model hence [0])
    AFprimarySeq = get1LetterPrimarySequenceFromModel(AF_PDB_obj[0])
    
    # Check for perfect match
    if AFprimarySeq == sourceSequence:
        if debug:
            print("Perfect Match")
        return {'NoMatch': False,
                'ExactMatch': True, 
                'TruncatedMatch': False,
                'AFTruncated':False,
                'SourceTruncated':False,
                'AFStartResidueOverlap': 1, 
                'AFEndResidueOverlap':len(AFprimarySeq),
                'SourceStartResidueOverlap':1, 
                'SourceEndResidueOverlap':len(sourceSequence)}

    # AF truncated sequence, but perfect match for what exists
    elif len(AFprimarySeq) < len(sourceSequence) and AFprimarySeq in sourceSequence:
        ResidueNumberOfMatch = sourceSequence.find(AFprimarySeq) + 1
        if debug:
            print(f'AF Truncated relative to Source: AF {len(AFprimarySeq)} residue sequence matched in Source\'s {len(sourceSequence)} residue sequence, starting at res number {ResidueNumberOfMatch} in Source')
        return {'NoMatch': False,
                'ExactMatch': False, 
                'TruncatedMatch': True,
                'AFTruncated':True,
                'SourceTruncated':False,
                'AFStartResidueOverlap': 1, 
                'AFEndResidueOverlap':len(AFprimarySeq),
                'SourceStartResidueOverlap':1, 
                'SourceEndResidueOverlap':len(AFprimarySeq)}
    
    # Target truncated sequnce, but perfect match for what exists
    elif len(AFprimarySeq) > len(sourceSequence) and sourceSequence in AFprimarySeq:
        ResidueNumberOfMatch = AFprimarySeq.find(sourceSequence) + 1
        if debug:
            print(f'Source sequence Truncated relative to AF: Source {len(sourceSequence)} residue sequence matched in AF\'s {len(AFprimarySeq)} residue sequence, starting at res number {ResidueNumberOfMatch} in AF')
        return {'NoMatch': False,
                'ExactMatch': False, 
                'TruncatedMatch': True,
                'AFTruncated':False,
                'SourceTruncated':True,
                'AFStartResidueOverlap': 1, 
                'AFEndResidueOverlap':len(sourceSequence),
                'SourceStartResidueOverlap':1, 
                'SourceEndResidueOverlap':len(sourceSequence)}
    # Fuzzy Matching - TODO
    # No Match
    else:
        if debug:
            print('NO MATCH')
        return utility_getNoSequenceMathcDict()

def checkAFPDBSequenceForDataFrame(dataframe, pathToPDBRoot, pdbParser, pdbPathFeatureName='PDB_path', sequenceFeatureName='Sequence', debug=False):
    """
    Reads dataframe containing sequences and paths to pdbs (generated from FetchAlphaFoldPDBs)
    Checks if each source sequence for AlphaFold match with  checkAFPDBSequence    
    Paths in AF info csvs are relatie so must also specifiy a path to PDB root dir
    returns dataframe with information on sequence match
    """
    matchInfo = []
    for index, row in dataframe.iterrows():
        sourceSequence = row[sequenceFeatureName]
        relative_pdb_path = row[pdbPathFeatureName]
        
        if not isinstance(relative_pdb_path, str):
            # Catch where no path exists
            _matchDict = utility_getNoPDBDict()
        else: 
            # Get full path (paths in csv are relative)
            fullPDBPath = os.path.join(pathToPDBRoot, row['PDB_path'])

            # Read in PDB
            structure = pdbParser.get_structure(str(index), fullPDBPath)

            # Check Match
            _matchDict = checkAFPDBSequence(sourceSequence, structure, debug=debug)
        
        _matchDict[sequenceFeatureName] = sourceSequence
        matchInfo.append(_matchDict)

    # Join Match Info to Original DF
    return pd.merge(dataframe, pd.DataFrame(matchInfo), left_on=sequenceFeatureName, right_on=sequenceFeatureName)

# ================================================================================================================================
#   RUN (End to end)
# ================================================================================================================================

def CheckAlphaFoldPDBSequences_EndToEnd(CONFIG):
    # Unpack CONFIG file
    SOURCE_SEQUENCES_PATH, ID_FEATURE_NAME, SEQUENCE_FEATURE_NAME, AF_PDB_DATAFRAME, AF_PDB_ROOT_DIR, OUTPUT_SAVEPATH = utility_UnPackConfig(CONFIG)

    # Get PDB parser
    pdbParser = PDBParser()

    # Load Source and PDM paths
    sourceProteinData = pd.read_csv(SOURCE_SEQUENCES_PATH)
    afPDBData = pd.read_csv(AF_PDB_DATAFRAME)
    print(f'{len(sourceProteinData)} protein sequences loaded, {len(afPDBData)} AF PDB entries to check')

    # Merge
    merged = sourceProteinData.merge(afPDBData, left_on=ID_FEATURE_NAME, right_on='uniprot_ID_source')
    assert len(merged)==len(sourceProteinData), "Merging with PDB information has led to duplications! Likely PDB info csv contains too many uniprot IDs"

    # Check Sequences
    print('Checking Sequences')
    matchedDF = checkAFPDBSequenceForDataFrame(merged, AF_PDB_ROOT_DIR, pdbParser, sequenceFeatureName=SEQUENCE_FEATURE_NAME)

    # Save output
    matchedDF.to_csv(OUTPUT_SAVEPATH, index=False)

    # Small debug output
    numSequences = len(matchedDF)
    numExactMatches = sum(matchedDF['ExactMatch'].values)
    numTruncatedMatches = sum(matchedDF['TruncatedMatch'].values)
    
    PDBs = matchedDF['PDB_path']
    noPDBs = [m for m in PDBs if not isinstance(m, str)]
    PDB_mask = [isinstance(m, str) for m in matchedDF['PDB_path']]
    numNoMatchInPDB = sum(matchedDF[PDB_mask]['NoMatch'].values)
    print(f'{numSequences} Sequences ran\n\t{numExactMatches} exact sequence matches in PDB\n\t{numTruncatedMatches} truncated sequence matches in PDB\n\t{numNoMatchInPDB} no sequence match in PDB\n\t{len(noPDBs)} PDB available')


# ================================================================================================================================
#   UTILITIES (Small functions either to tidy code for readability, or of limited/specific use)
# ================================================================================================================================

def utility_getNoPDBDict():
    return {'NoMatch': True,
            'ExactMatch': False, 
            'TruncatedMatch': False,
            'AFTruncated':False,
            'SourceTruncated':False,
            'AFStartResidueOverlap': 0, 
            'AFEndResidueOverlap': 0,
            'SourceStartResidueOverlap': 0, 
            'SourceEndResidueOverlap': 0}

def utility_getNoSequenceMathcDict():
    return utility_getNoPDBDict()

def utility_UnPackConfig(CONFIG):
    SOURCE_SEQUENCES_PATH = CONFIG['SOURCE_SEQUENCES_PATH']
    ID_FEATURE_NAME = CONFIG['ID_FEATURE_NAME']
    SEQUENCE_FEATURE_NAME = CONFIG['SEQUENCE_FEATURE_NAME']
    AF_PDB_DATAFRAME = CONFIG['AF_PDB_DATAFRAME']
    AF_PDB_ROOT_DIR = CONFIG['AF_PDB_ROOT_DIR']
    OUTPUT_SAVEPATH = CONFIG['OUTPUT_SAVEPATH']
    return SOURCE_SEQUENCES_PATH, ID_FEATURE_NAME, SEQUENCE_FEATURE_NAME, AF_PDB_DATAFRAME, AF_PDB_ROOT_DIR, OUTPUT_SAVEPATH

