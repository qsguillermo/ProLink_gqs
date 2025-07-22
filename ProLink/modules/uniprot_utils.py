import requests

def extract_protein_name(protein_data):
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
    return None

def format_protein_name_for_matching(protein_name):
    """
    Convierte el nombre de proteína para que coincida mejor con los labels:
    - Espacios → guiones bajos
    - Palabras en Title Case → minúsculas
    """
    if not protein_name:
        return ""
    parts = protein_name.split()
    formatted_parts = [
        part.lower() if part.istitle() else part
        for part in parts
    ]
    return "_".join(formatted_parts)

def get_protein_name_from_wp(wp_code):
    """
    Dado un código WP, devuelve el nombre formateado de la proteína.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "fields": "accession,protein_name",
        "query": f"xref:RefSeq-{wp_code}",
        "format": "json"
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        results = data.get("results", [])
        if not results:
            print(f"⚠️ No se encontró ninguna entrada para {wp_code}")
            return ""
        protein_data = results[0].get("proteinDescription", {})
        raw_name = extract_protein_name(protein_data)
        return format_protein_name_for_matching(raw_name)
    except Exception as e:
        print(f"❌ Error al consultar UniProt para {wp_code}: {e}")
        return ""
