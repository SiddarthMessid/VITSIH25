"""
PDB structure utilities for downloading and processing structures
"""

import os
import urllib.request
import urllib.error
from Bio.PDB.PDBParser import PDBParser
import requests


def download_pdb_structure(pdb_id, download_dir="pdb_files"):
    """
    Download PDB structure from RCSB PDB
    
    Args:
        pdb_id: PDB identifier (4-character code)
        download_dir: Directory to save PDB files (default: "pdb_files")
    
    Returns:
        File path if successful, None if failed
    """
    # Create directory if it doesn't exist
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    
    pdb_id = pdb_id.upper()
    file_path = os.path.join(download_dir, f"{pdb_id}.pdb")
    
    # Check if file already exists
    if os.path.exists(file_path):
        return file_path
    
    try:
        # Try multiple sources for PDB files
        urls = [
            f"https://files.rcsb.org/download/{pdb_id}.pdb",
            f"https://www.rcsb.org/pdb/files/{pdb_id}.pdb"
        ]
        
        for url in urls:
            try:
                urllib.request.urlretrieve(url, file_path)
                # Verify the file was downloaded and is valid
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    return file_path
            except (urllib.error.HTTPError, urllib.error.URLError):
                continue
        
        # If all URLs failed, try requests
        response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        if response.status_code == 200:
            with open(file_path, 'w') as f:
                f.write(response.text)
            return file_path
        
        return None
        
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None


def get_structure_info(pdb_file_path):
    """
    Extract basic information from PDB structure
    
    Args:
        pdb_file_path: Path to PDB file
    
    Returns:
        Dictionary containing structure information
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file_path)
        
        if structure is None:
            return {}
        
        # Get basic structure information
        info = {
            'num_models': len(list(structure.get_models())),
            'num_chains': len(list(structure.get_chains())),
            'num_residues': len(list(structure.get_residues())),
            'num_atoms': len(list(structure.get_atoms())),
            'chains': [chain.id for chain in structure.get_chains()]
        }
        
        return info
        
    except Exception as e:
        print(f"Error parsing PDB file {pdb_file_path}: {e}")
        return {}


def get_pdb_sequence(pdb_file_path, chain_id='A'):
    """
    Extract protein sequence from PDB file
    
    Args:
        pdb_file_path: Path to PDB file
        chain_id: Chain identifier (default: 'A')
    
    Returns:
        Protein sequence as string
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file_path)
        
        # Standard amino acid mapping
        aa_dict = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        
        sequence = ""
        
        if structure is not None:
            for model in structure:
                for chain in model:
                    if chain.id == chain_id:
                        for residue in chain:
                            if residue.id[0] == ' ':  # Standard residue
                                resname = residue.get_resname()
                                if resname in aa_dict:
                                    sequence += aa_dict[resname]
                                else:
                                    sequence += 'X'  # Unknown amino acid
                        break
                break
        
        return sequence
        
    except Exception as e:
        print(f"Error extracting sequence from {pdb_file_path}: {e}")
        return ""


def validate_pdb_file(pdb_file_path):
    """
    Validate PDB file format and content
    
    Args:
        pdb_file_path: Path to PDB file
    
    Returns:
        Boolean indicating if file is valid
    """
    try:
        if not os.path.exists(pdb_file_path):
            return False
        
        if os.path.getsize(pdb_file_path) == 0:
            return False
        
        # Check if file contains PDB header
        with open(pdb_file_path, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        # Look for common PDB record types
        pdb_records = ['HEADER', 'ATOM', 'HETATM', 'REMARK', 'SEQRES']
        
        for line in first_lines:
            if any(line.startswith(record) for record in pdb_records):
                return True
        
        return False
        
    except Exception:
        return False
