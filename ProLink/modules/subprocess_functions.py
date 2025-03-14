
import logging
import subprocess
import re

from .. import ProLink_path

logger = logging.getLogger()

def clean_newick_labels(newick_str, protein_name='alkene_reductase'):
    """
    Cleans the labels in a Newick tree by removing:
      - WP codes (pattern: WP_\d{9}\.\d)
      - The word "MULTISPECIES:" if present
      - The protein name (default 'alkene_reductase')
      - The trailing "Same_Domains" (and similar variants with hyphens or spaces)
    
    Returns the cleaned Newick string.
    """
    # Remove WP codes (e.g., WP_058328214.1)
    newick_str = re.sub(r'WP_\d{9}\.\d', '', newick_str)
    # Remove "MULTISPECIES:" if present
    newick_str = re.sub(r'MULTISPECIES:\s*', '', newick_str, flags=re.IGNORECASE)
    # Remove the protein name
    newick_str = re.sub(protein_name, '', newick_str, flags=re.IGNORECASE)
    # Remove trailing "Same_Domains" (with potential hyphens or spaces)
    newick_str = re.sub(r'[-_]*Same[_\s]*Domains', '', newick_str, flags=re.IGNORECASE)
    # Replace multiple underscores with a single underscore and clean extra spaces
    newick_str = re.sub(r'_+', '_', newick_str)
    newick_str = re.sub(r'\s+', ' ', newick_str)
    return newick_str.strip()


def align(muscle_input:str, muscle_output:str) -> None:
    '''
    Run a local alignment with MUSCLE v5

    Parameters
    ----------
    muscle_input : str
        Path of the input MUSCLE file
    muscle_output : str
        Path of the output MUSCLE file
    '''
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str) -> None:
    '''
    Run MEGA-CC to generate a phylogenetic tree

    Parameters
    ----------
    tree_type : str
        Type of tree to generate
    bootstrap_replications : int
        Number of bootstrap replications
    muscle_output : str
        Path of the input file (FASTA format from MUSCLE)
    mega_output : str
        Path of the MEGA-CC output file
    '''
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")
    mega_run = subprocess.run(mega_cmd)
    if mega_run.returncode != 0:
        logger.error(f"ERROR: MEGA-CC failed")
        raise RuntimeError(f"MEGA-CC failed")

    # Read the generated Newick tree and clean its labels
    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        cleaned_newick = clean_newick_labels(newick)
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise
