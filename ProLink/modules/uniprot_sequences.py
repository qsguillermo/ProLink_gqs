
import requests
import re
import logging
from Bio import SeqIO  # Para manejar archivos FASTA correctamente

logger = logging.getLogger()

def check_uniprot_batch(wp_codes):
    """
    Verifica la existencia de múltiples códigos WP en UniProt en una sola consulta.

    Parámetros:
    wp_codes (list): Lista de códigos WP a verificar.

    Retorna:
    set: Conjunto de códigos WP encontrados en UniProt.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    query = " OR ".join(wp_codes)  # Unimos todos los códigos en una sola consulta
    params = {"query": query, "fields": "accession", "format": "json"}

    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            return {entry['primaryAccession'] for entry in data.get("results", [])}
        else:
            logger.warning(f"Error en la consulta a UniProt: {response.status_code}")
            return set()
    except Exception as e:
        logger.error(f"Error al conectar con UniProt: {e}")
        return set()

def filter_valid_sequences(input_fasta, output_fasta):
    """
    Filtra las secuencias eliminando aquellas cuyos códigos WP no existen en UniProt.
    Si una secuencia no tiene código WP_, se conserva.

    Parámetros:
    input_fasta (str): Archivo FASTA de entrada con secuencias.
    output_fasta (str): Archivo FASTA de salida con secuencias válidas.
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Extraer códigos WP_ de las descripciones
    wp_data = {}
    for seq in sequences:
        match = re.search(r'(WP_\d{9}\.\d)', seq.description)
        if match:
            wp_data[seq.id] = match.group(1)

    logger.info(f"Número total de secuencias: {len(sequences)}")
    logger.info(f"Número de códigos WP encontrados: {len(wp_data)}")

    # Verificar en UniProt en lotes de 100
    wp_codes = list(set(wp_data.values()))  # Evitar códigos duplicados
    valid_wp_codes = set()
    batch_size = 100  

    for i in range(0, len(wp_codes), batch_size):
        batch = wp_codes[i:i+batch_size]
        valid_wp_codes.update(check_uniprot_batch(batch))

    # Filtrar secuencias válidas
    valid_sequences = [
        seq for seq in sequences 
        if (seq.id in wp_data and wp_data[seq.id] in valid_wp_codes) or seq.id not in wp_data
    ]

    # Guardar el nuevo archivo FASTA solo con secuencias válidas
    SeqIO.write(valid_sequences, output_fasta, "fasta")

    logger.info(f"Secuencias válidas después del filtrado: {len(valid_sequences)}")
    logger.info(f"Resultados guardados en {output_fasta}")
