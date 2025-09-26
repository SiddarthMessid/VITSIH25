"""
BLAST search utilities for protein sequence analysis
"""

import time
from Bio.Blast import NCBIWWW
from Bio import SearchIO
from io import StringIO


def perform_blast_search(sequence, program="blastp", database="pdb", expect=0.001, hitlist_size=20):
    """
    Perform BLAST search against PDB database
    
    Args:
        sequence: Protein sequence to search
        program: BLAST program (default: blastp)
        database: Database to search (default: pdb)
        expect: E-value threshold (default: 0.001)
        hitlist_size: Maximum number of hits to return (default: 20)
    
    Returns:
        BLAST results object or None if failed
    """
    try:
        # Perform BLAST search
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=sequence,
            expect=expect,
            hitlist_size=hitlist_size
        )
        
        # Parse BLAST results
        blast_results = SearchIO.read(result_handle, "blast-xml")
        result_handle.close()
        
        return blast_results
        
    except Exception as e:
        print(f"Error performing BLAST search: {e}")
        return None


def extract_pdb_ids_from_blast(blast_results, top_n=10):
    """
    Extract PDB IDs and relevant information from BLAST results
    
    Args:
        blast_results: BLAST results object from SearchIO
        top_n: Number of top hits to extract (default: 10)
    
    Returns:
        List of dictionaries containing PDB information
    """
    pdb_list = []
    
    for i, hit in enumerate(blast_results[:top_n]):
        # Extract PDB ID from hit ID
        hit_id = hit.id
        
        # Handle different PDB ID formats
        if '|' in hit_id:
            # Format: pdb|7D4F|A or similar
            parts = hit_id.split('|')
            pdb_id = parts[1][:4] if len(parts) > 1 else hit_id[:4]
            chain = parts[2] if len(parts) > 2 else 'A'
        else:
            # Direct PDB ID format
            pdb_id = hit_id[:4]
            chain = 'A'
        
        # Get best HSP (High-scoring Segment Pair)
        best_hsp = hit[0]
        
        pdb_info = {
            'rank': i + 1,
            'pdb_id': pdb_id.upper(),
            'chain': chain,
            'full_id': hit_id,
            'description': hit.description,
            'evalue': best_hsp.evalue,
            'bitscore': best_hsp.bitscore,
            'identity': best_hsp.ident_num,
            'alignment_length': best_hsp.aln_span,
            'query_start': best_hsp.query_start,
            'query_end': best_hsp.query_end,
            'hit_start': best_hsp.hit_start,
            'hit_end': best_hsp.hit_end,
            'query_seq': str(best_hsp.query.seq),
            'hit_seq': str(best_hsp.hit.seq),
            'match_seq': best_hsp.aln_annotation.get('similarity', '')
        }
        
        pdb_list.append(pdb_info)
    
    return pdb_list


def get_blast_statistics(pdb_hits):
    """
    Calculate statistics from BLAST results
    
    Args:
        pdb_hits: List of PDB hit dictionaries
    
    Returns:
        Dictionary containing statistics
    """
    if not pdb_hits:
        return {}
    
    evalues = [hit['evalue'] for hit in pdb_hits]
    bitscores = [hit['bitscore'] for hit in pdb_hits]
    identities = [hit['identity'] for hit in pdb_hits]
    
    stats = {
        'total_hits': len(pdb_hits),
        'min_evalue': min(evalues),
        'max_evalue': max(evalues),
        'avg_evalue': sum(evalues) / len(evalues),
        'min_bitscore': min(bitscores),
        'max_bitscore': max(bitscores),
        'avg_bitscore': sum(bitscores) / len(bitscores),
        'min_identity': min(identities),
        'max_identity': max(identities),
        'avg_identity': sum(identities) / len(identities)
    }
    
    return stats
