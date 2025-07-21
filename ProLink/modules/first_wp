import requests

def get_wp_from_code(code: str) -> str:
    """
    Obtiene el código WP a partir de un código de entrada (como EMBL).
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"({code})",
        "fields": "accession",
        "format": "json"
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        if not data.get("results"):
            return None

        accession = data["results"][0].get("primaryAccession", None)
        if not accession:
            return None

        # Obtener el WP desde la entrada
        entry_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        entry_resp = requests.get(entry_url)
        entry_resp.raise_for_status()
        entry = entry_resp.json()

        for ref in entry.get("uniProtKBCrossReferences", []):
            if ref.get("database") == "RefSeq":
                wp = ref.get("id")
                if wp and wp.startswith("WP_"):
                    return wp
        return None

    except Exception as e:
        print(f"❌ Error al obtener WP desde {code}: {e}")
        return None
