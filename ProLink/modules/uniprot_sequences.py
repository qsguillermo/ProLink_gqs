import requests
import re
import logging
from Bio import SeqIO  # To properly handle FASTA files

logger = logging.getLogger()

url = "https://rest.uniprot.org/uniprotkb/search"

def check_uniprot_single(wp_code):
    """
    Verify the existence of a single WP code in UniProt.
    
    Parameters:
    wp_code (str): WP code to verify.

    Returns:
    bool: True if the WP code exists in UniProt, False otherwise.
    """
    params = {
        "query": f"xref:RefSeq-{wp_code}",
        "fields": "accession",
        "format": "json",
        "size": 1  # We only need to check if it exists
    }

    print(f"Consulta a UniProt: {params['query']}")  # Debug: Check what is being sent to UniProt
    
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        return bool(data.get("results"))  # Returns True if there's at least one result
    except requests.exceptions.RequestException as e:
        logger.error(f"Error al conectar con UniProt: {e}")
        return False

def filter_valid_sequences(input_fasta, output_fasta):
    """
    Filters sequences by removing those whose WP codes do not exist in UniProt.
    Sequences without a WP_ code are retained.

    Parameters:
    input_fasta (str): Input FASTA file with sequences.
    output_fasta (str): Output FASTA file with valid sequences.
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Extract WP_ codes from sequence descriptions
    wp_data = {}
    for seq in sequences:
        match = re.search(r'((?:WP|XP|NP)_\d{9}\.\d)', seq.description)
        if match:
            wp_data[seq.description] = match.group(1)
    
    print(f"Códigos WP extraídos: {list(wp_data.values())}")  # Debug: Show extracted WP codes
    
    logger.info(f"Número total de secuencias: {len(sequences)}")
    logger.info(f"Número de códigos WP encontrados: {len(wp_data)}")

    # Verify each WP code in UniProt individually
    valid_wp_codes = {wp for wp in wp_data.values() if check_uniprot_single(wp)}

    # Filter valid sequences
    valid_sequences = [
        seq for seq in sequences 
        if seq.description not in wp_data or wp_data[seq.description] in valid_wp_codes
    ]
   
    # Write the valid sequences to the new FASTA file
    SeqIO.write(valid_sequences, output_fasta, "fasta")
    print(f"Secuencias válidas después del filtrado: {len(valid_sequences)}")  # Debug: Show number of valid sequences
    logger.info(f"Resultados guardados en {output_fasta}")

import csv

wp_codes = list(wp_data)

output_file = "annotation.csv"
results = []

for wp in wp_codes:
    query_string = f"({wp})"  # búsqueda libre
    params = {
        "fields": "accession",
        "query": query_string,
        "format": "json"
    }

    response = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params)
    data = response.json()

    if data.get("results"):
        for r in data["results"]:
            primary_accession = r.get("primaryAccession", "Not found")
            results.append({
                "WP_code": wp,
                "UniProt_primaryAccession": primary_accession
            })
    else:
        results.append({
            "WP_code": wp,
            "UniProt_primaryAccession": "Not found"
        })

# Escribir resultados al CSV
with open(output_file, mode="w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=["WP_code", "UniProt_primaryAccession"], delimiter=';')
    writer.writeheader()
    writer.writerows(results)
    logger.info (f"\n✅ Archivo CSV_annotation generado: {output_file}")

