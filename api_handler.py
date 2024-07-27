import requests

class AlphaFold_DB_API():
    """
    A class for interacting with the AlphaFold Database API.

    Attributes:
        uniprot_accession (str): The UniProt accession code of the protein.
        result (dict or None): The JSON result from the API, or None if unavailable.

    Methods:
        get: Retrieves data from the AlphaFold Database API.
    """
    def __init__(self, uniprot_accession):
        self.uniprot_accession = uniprot_accession
        self.result = None  # Initialize result to None

    def get(self):
        """Retrieves data from the AlphaFold Database API.

        Returns:
            dict or None: The JSON result from the API, or None if an error occurs.
        """
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{self.uniprot_accession}"
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for HTTP errors
            self.result = response.json()
            return self.result
        except requests.exceptions.RequestException as e:
            print(f"Error retrieving data for {self.uniprot_accession}: {e}")
            return None

class retrieve_AFDB(AlphaFold_DB_API):
    """
    A class for retrieving AlphaMissense and PDB URLs from the AlphaFold Database.

    Attributes:
        uniprot_accession (str): The UniProt accession code of the protein.
        result (dict or None): The JSON result from the API, or None if unavailable.

    Methods:
        get_AM_url: Retrieves the AlphaMissense URL, if available.
        get_pdb_url: Retrieves the PDB URL, if available.
        get_PAE_url: Retrieves the PAE image.
    """
    def __init__(self, uniprot_accession):
        super().__init__(uniprot_accession)  # Call the parent constructor
        self.get()  # Attempt to retrieve the data

    def get_AM_url(self):
        """Retrieves the AlphaMissense URL, if available.

        Returns:
            str: The AlphaMissense URL, or an error message if unavailable.
        """
        if self.result:
            return self.result[0].get('amAnnotationsUrl',
                    f"The UniProt ID {self.uniprot_accession} does not contain data for AlphaMissense.")
        else:
            return f"Error retrieving data for {self.uniprot_accession}"

    def get_pdb_url(self):
        """Retrieves the PDB URL, if available.

        Returns:
            str: The PDB URL, or an error message if unavailable.
        """
        if self.result:
            return self.result[0].get('pdbUrl', "Failed to retrieve PDB file from AlphaFold API.")
        else:
            return f"Error retrieving data for {self.uniprot_accession}"


    def get_pae_url(self):
            """Retrieves the PAE URL, if available.

            Returns:
                str: The PAE URL, or an error message if unavailable.
            """
            if self.result:
                return self.result[0].get('paeImageUrl', "Failed to retrieve PAE from AlphaFold API.")
            else:
                return f"Error retrieving data for {self.uniprot_accession}"
