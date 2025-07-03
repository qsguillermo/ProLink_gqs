import requests
import re
import logging
import csv
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
    logger.info(f"Número de códigos WP_ encontrados: {len(wp_data)}")

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

    try:
        annotate_uniprot_codes(valid_wp_codes)  # Call the annotation function
        print("annotate_uniprot_codes completed successfully.")
    except Exception as e:
        logger.warning(f"annotate_uniprot_codes failed: {e}")
        print(f"Warning: Anotacion fallida: {e}")

ef extract_protein_name(protein_data):
    """
    Extrae el mejor nombre posible para la proteína:
    1. recommendedName
    2. submissionNames[0]
    3. alternativeNames[0]
    """
    try:
        return protein_data["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        pass
    try:
        return protein_data["submissionNames"][0]["fullName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    try:
        return protein_data["alternativeNames"][0]["fullName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    return "Not found"

def extract_ec_number(protein_data):
    """
    Extrae el primer número EC si está presente.
    """
    try:
        ec_list = protein_data.get("recommendedName", {}).get("ecNumbers", [])
        if ec_list:
            return ec_list[0]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    return "None"

def get_cofactors_from_accession(accession):
    """
    Consulta UniProt con accession para extraer cofactores dentro de comments.
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
    except Exception as e:
        print(f"❌ Error al obtener entrada UniProt para {accession}: {e}")
        return "Error"

    cofactors = []

    comments = entry.get("comments")
    if comments is None:
        print(f"⚠️ 'comments' no existe para {accession}")
        return "None"
    if not isinstance(comments, list):
        print(f"⚠️ 'comments' no es una lista para {accession}. Es: {type(comments)}")
        return "None"

    for comment in comments:
        if comment.get("commentType") == "COFACTOR":
            for cofactor in comment.get("cofactors", []):
                name_field = cofactor.get("name")
                if isinstance(name_field, dict):
                    name = name_field.get("value")
                elif isinstance(name_field, str):
                    name = name_field
                else:
                    name = None
                if name:
                    cofactors.append(name)


    return "; ".join(cofactors) if cofactors else "None"

def annotate_uniprot_codes(valid_wp_codes, output_file="annotationE4.csv"):
    results = []

    for wp in valid_wp_codes:
        query_string = f"xref:RefSeq-{wp}"
        params = {
            "fields": "accession,organism_name,protein_name",
            "query": query_string,
            "format": "json"
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()

            if data.get("results"):
                for r in data["results"]:
                    print(f"DEBUG r: {r}")  # DEBUG

                    protein_data = r.get("proteinDescription")
                    if not isinstance(protein_data, dict):
                        protein_data = {}

                    protein_name = extract_protein_name(protein_data)
                    ec_number = extract_ec_number(protein_data)
                    accession = r.get("primaryAccession", "Not found")

                    cofactors = get_cofactors_from_accession(accession) if accession != "Not found" else "None"

                    results.append({
                        "WP_code": wp,
                        "UniProt_accession": accession,
                        "Organism": r.get("organism", {}).get("scientificName", "Not found"),
                        "Protein_name": protein_name,
                        "EC_number": ec_number,
                        "Cofactors": cofactors
                    })
            else:
                results.append({
                    "WP_code": wp,
                    "UniProt_accession": "Not found",
                    "Organism": "Not found",
                    "Protein_name": "Not found",
                    "EC_number": "None",
                    "Cofactors": "None"
                })
        except Exception as e:
            print(f"❌ Error al consultar {wp}: {e}")
            results.append({
                "WP_code": wp,
                "UniProt_accession": "error",
                "Organism": "error",
                "Protein_name": "error",
                "EC_number": "error",
                "Cofactors": "error"
            })

    with open(output_file, mode="w", newline="", encoding="utf-8") as file:
        fieldnames = ["WP_code", "UniProt_accession", "Organism", "Protein_name", "EC_number", "Cofactors"]
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()
        writer.writerows(results)

    print(f"✅ Archivo CSV generado: {output_file}")
