
import logging
import subprocess
import time
import os
import re

from .. import ProLink_path

logger = logging.getLogger()

def clean_label(label, protein_name='alkene_reductase'):
    """
    Cleans a single Newick label by removing:
      - WP codes (pattern: WP_\d{9}\.\d)
      - The term "MULTISPECIES:" if present
      - The protein name (default 'alkene_reductase')
      - Trailing "Same Domains" (or similar variants)
    Then, it extracts and returns only the species name and the cluster marker,
    abbreviating the genus name if possible.
    """
    label = label.strip("'\"")
    label = re.sub(r"WP[\s_]\d{9}\.\d", "", label)
    label = re.sub(r"MULTISPECIES:\s*", "", label, flags=re.IGNORECASE)
    protein_regex = re.escape(protein_name).replace(r'\_', r'[\s_]+')
    label = re.sub(protein_regex, "", label, flags=re.IGNORECASE)
    label = re.sub(r"\bunclassified\b", "", label, flags=re.IGNORECASE)
    label = re.sub(r"[-_]*Same[_\s]*Domains", "", label, flags=re.IGNORECASE)
    label = label.strip(" _")

    pattern = re.compile(
        r"(?P<species>[A-Z][a-zA-Z0-9_\/]*(?:[\s_]+[a-zA-Z0-9_\/\.\-]+)*)"
        r"[\s_-]+(?P<cluster>---C\d+)",
        flags=re.IGNORECASE
    )
    match = pattern.search(label)
    if match:
        species = match.group('species').strip().replace(" ", "_")
        cluster = match.group('cluster').strip()

        # Abreviar el género si es posible y no es "sp."
        parts = species.split("_")
        if len(parts) >= 2 and not parts[1].startswith("sp"):
            genus_initial = parts[0][0]
            species_abbr = f"{genus_initial}_{'_'.join(parts[1:])}"
        else:
            species_abbr = species

        return f"{species_abbr}_{cluster}"
    else:
        return label

def clean_newick_string(newick_str, protein_name='alkene_reductase'):
    pattern = re.compile(
        r"('([^']+---C\d+[^']*)'|\"([^\"]+---C\d+[^\"]*)\"|([A-Za-z0-9 _\.\-]+---C\d+))",
        flags=re.IGNORECASE
    )
    def replacer(match):
        full_label = match.group(0)
        return clean_label(full_label, protein_name)
    return pattern.sub(replacer, newick_str)

def align(muscle_input:str, muscle_output:str) -> None:
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str) -> None:
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")

    mega_run = subprocess.run(mega_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    logger.debug(f"MEGA-CC stdout: {mega_run.stdout}")
    logger.debug(f"MEGA-CC stderr: {mega_run.stderr}")
    if mega_run.returncode != 0:
        logger.error("ERROR: MEGA-CC failed")
        raise RuntimeError("MEGA-CC failed")

    time.sleep(5)

    if not os.path.exists(mega_output):
        alternative = mega_output.rsplit('.', 1)[0] + ".mega"
        if os.path.exists(alternative):
            logging.info(f"Using alternative output file: {alternative}")
            mega_output = alternative
        else:
            logger.error(f"ERROR: MEGA-CC did not produce the output file: {mega_output}")
            raise FileNotFoundError(f"Output file {mega_output} not found")

    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        cleaned_newick = clean_newick_string(newick, protein_name='alkene_reductase')
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
        logging.info("✅ Testing_abb")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise
