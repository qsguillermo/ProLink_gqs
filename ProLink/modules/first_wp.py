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
                wp_of_interest = ref.get("id")
                if wp_of_interest and wp_of_interest.startswith("WP_"):
                    return wp_of_interest
        return None

    except Exception as e:
        print(f"❌ Error al obtener WP desde {code}: {e}")
        return None

def reorder_fasta_with_study_sequence(txt_file, fasta_file, wp_of_interest, study_seq, output_file):
    cluster_dict = {}
    current_cluster = None

    # Paso 1: Leer archivo .txt de clústeres y construir diccionario
    with open(txt_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#Cluster"):
                current_cluster = "C" + line.split()[-1]
            elif line and current_cluster:
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    wp_code, descriptor = parts
                    full_header = f"{wp_code} {descriptor}"
                    cluster_dict[wp_code] = {
                        "cluster": current_cluster,
                        "full_header": full_header
                    }

    # Paso 2: Leer archivo FASTA original
    from collections import defaultdict

    fasta_entries = defaultdict(list)
    with open(fasta_file, "r") as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    wp = current_header.split()[0][1:]  # Extrae WP
                    fasta_entries[wp] = [current_header, "".join(current_seq)]
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and current_seq:
            wp = current_header.split()[0][1:]
            fasta_entries[wp] = [current_header, "".join(current_seq)]

    # Paso 3: Construir el nuevo FASTA reordenado
    clusters_written = set()
    with open(output_file, "w") as out:
        for wp, info in cluster_dict.items():
            cluster = info["cluster"]
            full_header = info["full_header"]

            if cluster not in clusters_written:
                # Si el WP de estudio pertenece a este clúster, escríbelo primero
                if cluster_dict.get(wp_of_interest, {}).get("cluster") == cluster:
                    study_header = f">{cluster_dict[wp_of_interest]['full_header']}---{cluster}"
                    out.write(study_header + "\n")
                    out.write(study_seq.strip() + "\n")
                
                clusters_written.add(cluster)

            if wp == wp_of_interest:
                continue  # Ya fue escrito antes

            if wp in fasta_entries:
                entry_header, seq = fasta_entries[wp]
                new_header = f">{full_header}---{cluster}"
                out.write(new_header + "\n")
                out.write(seq + "\n")

    print(f"✅ Archivo FASTA reordenado generado en: {output_file}")

txt_file = "seqs_cluster.txt"
fasta_file = "seqs_cluster.fasta"
wp_of_interest = "WP_006015261.1"
study_seq = "MSKMFTDLKRHGLNGPTMRTRWSAIFYAAADFDIAPLR..."  # tu secuencia aquí
output_file = "seqs_cluster_interest.fasta"

