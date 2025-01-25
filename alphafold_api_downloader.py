import requests
import logging
import argparse
import os
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def sleep_and_retry(func, timeout=5):
    def wrapper(*args, **kwargs):
        while True:
            try:
                return func(*args, **kwargs)
            except requests.exceptions.RequestException:
                time.sleep(timeout)

    return wrapper

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


def extract_urls(data, uniprot_accession, extract_am_files=False):
    """Extracts all files from the API data.

    Args:
        data (dict): The JSON data from the API.

    Returns:
        str: The AlphaMissense URL, or an error message.
    """
    if data:
        logging.info(f"Retrieving data for {uniprot_accession}")
        result = {}
        result['alphaMissenseUrl'] = data[0].get('amAnnotationsUrl', f"No AlphaMissense data for {uniprot_accession}")
        result['pdbUrl'] = data[0].get('pdbUrl', f"No PDB URL for {uniprot_accession}")
        result['cifUrl'] = data[0].get('cifUrl', f"No PDB URL for {uniprot_accession}") 
        result['paeImageUrl'] = data[0].get('paeImageUrl', f"No PDB URL for {uniprot_accession}") 
        
# Extract additional files if extract_am_files is True
        if extract_am_files:
            result['amAnnotationsHg19Url'] = data[0].get('amAnnotationsHg19Url', f"No amAnnotationsHg19Url data for {uniprot_accession}")
            result['amAnnotationsHg38Url'] = data[0].get('amAnnotationsHg38Url', f"No amAnnotationsHg38Url data for {uniprot_accession}")

        return result
    else:
        return {
            'alphaMissenseUrl': f"Error: No data provided for {uniprot_accession}",
            'pdbUrl': f"Error: No data provided for {uniprot_accession}"
        }


def download_file(url, output_dir):
    """
    Downloads a file from a given URL.

    Args:
        url (str): The URL of the file to download.
        output_dir (str): The directory to save the downloaded file.
    """
    
    if not url or not url.startswith("http"):  # Check for both None and invalid URLs
        return

    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()

        filename = os.path.join(output_dir, url.split('/')[-1])
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logging.info(f"Downloaded {filename}")

    except requests.exceptions.RequestException as e:
        logging.error(f"Error downloading {url}: {e}")



def process_uniprot_ids_from_file(filename, extract_am_files=False):
    """
    Reads UniProt IDs from a text file, fetches data and extracts URLs.

    Args:
        filename (str): The name of the file containing UniProt IDs.
        extract_am_files (bool): Whether to extract files from amAnnotationsHg19Url and amAnnotationsHg38Url.
    """
    
    with open(filename, 'r') as f:
        uniprot_ids = f.read().replace(' ', '').split(',')  # Remove whitespace and split by comma
        
    output_dir = "downloaded_files"
    os.makedirs(output_dir, exist_ok=True)

    for uniprot_id in uniprot_ids:
        try:
            data = sleep_and_retry(fetch_AFDB_data, timeout=5)(uniprot_id)
            urls = extract_urls(data, uniprot_id, extract_am_files)
            
            for url in urls.values():
                download_file(url, output_dir)
        
        except requests.exceptions.RequestException as e:
            print(f"Error processing {uniprot_id}: {e}")




# Define the main function
def main():
    parser = argparse.ArgumentParser(description='Fetch and extract URLs from the AlphaFold Database API.')
    parser.add_argument('filename', help='Text file containing UniProt IDs, one per line')
    parser.add_argument('--extract_am_files', action='store_true', help='Extract files from amAnnotationsHg19Url and amAnnotationsHg38Url (default: False)')  # Add argument for extract_am_files
    args = parser.parse_args()

    # Call the process_uniprot_ids_from_file function with the provided filename and output option
    process_uniprot_ids_from_file(args.filename, args.extract_am_files)

# Call the main function if this script is run directly
if __name__ == '__main__':
    main()