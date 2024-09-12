import requests
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def fetch_AFDB_data(uniprot_accession):
    """Fetches data from the AlphaFold Database API.

    Args:
        uniprot_accession (str): The UniProt accession code.

    Returns:
        dict or None: The JSON result, or None on error.
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
    try:
        logging.info(f"Fetching data for {uniprot_accession}")
        response = requests.get(url)
        response.raise_for_status()
        logging.info(f"Successfully retrieved data for {uniprot_accession}")
        return response.json()

    except requests.exceptions.RequestException as e:
        logging.error(f"Error retrieving data for {uniprot_accession}, check UniProt accession")
        return None


def extract_alpha_missense_url(data):
    """Extracts the AlphaMissense URL from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The AlphaMissense URL, or an error message.
    """
    if data:
        logging.info(f"Retrieving AM data")
        return data[0].get('amAnnotationsUrl', f"No AlphaMissense data for this protein.")
    else:
        return "Error: No data provided."

def extract_pdb_url(data):
    """Extracts the PDB URL from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The PDB url, or an error message.
    """
    if data:
        logging.info(f"Retrieving url for PDB file")
        return data[0].get('pdbUrl', "Failed to retrieve PDB URL.")
    else:
        return "Error: No data provided."

def extract_pae_url(data):
    """Extracts the PAE image URL from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The PAE image URL, or an error message.
    """
    if data:
        logging.info(f"Retrieving url for PAE png")
        return data[0].get('paeImageUrl', "Failed to retrieve PAE image URL.")
    else:
        return "Error: No data provided."