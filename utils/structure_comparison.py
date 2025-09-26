"""
Structure comparison utilities for overlaying and analyzing multiple protein structures
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import py3Dmol
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer


def create_structure_comparison_viewer(structure_data_list, width=700, height=500):
    """
    Create 3D viewer for comparing multiple protein structures
    
    Args:
        structure_data_list: List of dictionaries with keys:
            - 'pdb_path': Path to PDB file
            - 'pdb_id': PDB identifier
            - 'color': Color for the structure
            - 'label': Display label
        width: Viewer width (default: 700)
        height: Viewer height (default: 500)
    
    Returns:
        HTML string for embedding the viewer
    """
    try:
        if not structure_data_list:
            return "<p>No structures provided for comparison</p>"
        
        # Create py3Dmol viewer
        viewer = py3Dmol.view(width=width, height=height)
        
        # Color palette for structures
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow']
        
        # Add each structure with different colors
        for i, struct_data in enumerate(structure_data_list):
            pdb_path = struct_data['pdb_path']
            if not os.path.exists(pdb_path):
                continue
                
            # Read PDB file content
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
            
            # Use provided color or default from palette
            color = struct_data.get('color', colors[i % len(colors)])
            
            # Add model to viewer
            viewer.addModel(pdb_content, 'pdb')
            
            # Set style for this model
            viewer.setStyle(
                {'model': i}, 
                {'cartoon': {'color': color, 'opacity': 0.8}}
            )
            
            # Add label
            label = struct_data.get('label', struct_data.get('pdb_id', f'Structure {i+1}'))
            viewer.addLabel(
                label, 
                {'position': {'x': 0, 'y': 0, 'z': 0}, 'backgroundColor': color, 'fontSize': 12}
            )
        
        # Center and zoom to show all structures
        viewer.zoomTo()
        
        # Return HTML
        return viewer._make_html()
        
    except Exception as e:
        return f"<p>Error creating structure comparison viewer: {str(e)}</p>"


def calculate_structure_alignment(ref_pdb_path, target_pdb_path, chain_id='A'):
    """
    Calculate structural alignment between two protein structures
    
    Args:
        ref_pdb_path: Path to reference PDB file
        target_pdb_path: Path to target PDB file
        chain_id: Chain identifier to use for alignment
    
    Returns:
        Dictionary with alignment metrics (RMSD, aligned atoms count, etc.)
    """
    try:
        parser = PDBParser(QUIET=True)
        
        # Load structures
        ref_structure = parser.get_structure("reference", ref_pdb_path)
        target_structure = parser.get_structure("target", target_pdb_path)
        
        # Get chains
        ref_chain = None
        target_chain = None
        
        if ref_structure is not None:
            for model in ref_structure:
                for chain in model:
                    if chain.id == chain_id:
                        ref_chain = chain
                        break
                if ref_chain:
                    break
        
        if target_structure is not None:
            for model in target_structure:
                for chain in model:
                    if chain.id == chain_id:
                        target_chain = chain
                        break
                if target_chain:
                    break
        
        if not ref_chain or not target_chain:
            return {'error': f'Chain {chain_id} not found in one or both structures'}
        
        # Extract CA atoms for alignment
        ref_atoms = []
        target_atoms = []
        common_residues = []
        
        ref_residues = {res.id[1]: res for res in ref_chain if 'CA' in res}
        target_residues = {res.id[1]: res for res in target_chain if 'CA' in res}
        
        # Find common residue positions
        for res_id in sorted(set(ref_residues.keys()) & set(target_residues.keys())):
            if 'CA' in ref_residues[res_id] and 'CA' in target_residues[res_id]:
                ref_atoms.append(ref_residues[res_id]['CA'])
                target_atoms.append(target_residues[res_id]['CA'])
                common_residues.append(res_id)
        
        if len(ref_atoms) < 3:
            return {'error': 'Insufficient common CA atoms for alignment (need at least 3)'}
        
        # Perform structural superposition
        superimposer = Superimposer()
        superimposer.set_atoms(ref_atoms, target_atoms)
        
        # Calculate metrics
        rmsd = superimposer.rms
        aligned_atoms = len(ref_atoms)
        rotran = superimposer.rotran
        rotation_matrix = rotran[0] if rotran else None
        translation_vector = rotran[1] if rotran else None
        
        return {
            'rmsd': rmsd,
            'aligned_atoms': aligned_atoms,
            'common_residues': common_residues,
            'coverage_ref': len(common_residues) / len(ref_residues) * 100,
            'coverage_target': len(common_residues) / len(target_residues) * 100,
            'rotation_matrix': rotation_matrix.tolist() if rotation_matrix is not None else None,
            'translation_vector': translation_vector.tolist() if translation_vector is not None else None
        }
        
    except Exception as e:
        return {'error': f'Structure alignment failed: {str(e)}'}


def create_comparison_metrics_table(alignment_results):
    """
    Create a DataFrame with structure comparison metrics
    
    Args:
        alignment_results: List of alignment result dictionaries
    
    Returns:
        pandas DataFrame with comparison metrics
    """
    try:
        if not alignment_results:
            return pd.DataFrame()
        
        metrics_data = []
        
        for result in alignment_results:
            if 'error' not in result:
                metrics_data.append({
                    'Structure_Pair': result.get('pair_label', 'Unknown'),
                    'RMSD_Å': round(result['rmsd'], 3),
                    'Aligned_Atoms': result['aligned_atoms'],
                    'Coverage_Ref_%': round(result['coverage_ref'], 1),
                    'Coverage_Target_%': round(result['coverage_target'], 1),
                    'Quality': 'Excellent' if result['rmsd'] < 2.0 else 'Good' if result['rmsd'] < 5.0 else 'Poor'
                })
            else:
                metrics_data.append({
                    'Structure_Pair': result.get('pair_label', 'Unknown'),
                    'RMSD_Å': 'Error',
                    'Aligned_Atoms': 'N/A',
                    'Coverage_Ref_%': 'N/A',
                    'Coverage_Target_%': 'N/A',
                    'Quality': result['error']
                })
        
        return pd.DataFrame(metrics_data)
        
    except Exception as e:
        # Return empty DataFrame on error
        return pd.DataFrame({'Error': [f'Failed to create metrics table: {str(e)}']})


def create_rmsd_heatmap(rmsd_matrix, structure_labels):
    """
    Create a heatmap visualization of RMSD values between structures
    
    Args:
        rmsd_matrix: 2D numpy array with RMSD values
        structure_labels: List of structure labels
    
    Returns:
        Plotly figure object
    """
    try:
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=rmsd_matrix,
            x=structure_labels,
            y=structure_labels,
            colorscale='viridis_r',  # Reversed viridis (low RMSD = bright)
            hoverongaps=False,
            text=np.round(rmsd_matrix, 2),
            texttemplate="%{text} Å",
            textfont={"size": 10},
        ))
        
        fig.update_layout(
            title="Structural Similarity Matrix (RMSD)",
            xaxis_title="Structures",
            yaxis_title="Structures",
            width=600,
            height=500
        )
        
        return fig
        
    except Exception as e:
        # Return empty figure with error message
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating RMSD heatmap: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig


def perform_pairwise_alignment(structure_paths, structure_labels, chain_id='A'):
    """
    Perform pairwise structural alignments between all structures
    
    Args:
        structure_paths: List of PDB file paths
        structure_labels: List of structure labels
        chain_id: Chain identifier to use for alignment
    
    Returns:
        RMSD matrix and alignment results list
    """
    try:
        n_structures = len(structure_paths)
        rmsd_matrix = np.zeros((n_structures, n_structures))
        alignment_results = []
        
        for i in range(n_structures):
            for j in range(i+1, n_structures):
                result = calculate_structure_alignment(
                    structure_paths[i], 
                    structure_paths[j], 
                    chain_id
                )
                
                result['pair_label'] = f"{structure_labels[i]} vs {structure_labels[j]}"
                alignment_results.append(result)
                
                if 'error' not in result:
                    rmsd_matrix[i, j] = result['rmsd']
                    rmsd_matrix[j, i] = result['rmsd']
                else:
                    rmsd_matrix[i, j] = np.nan
                    rmsd_matrix[j, i] = np.nan
        
        return rmsd_matrix, alignment_results
        
    except Exception as e:
        return np.array([]), [{'error': f'Pairwise alignment failed: {str(e)}'}]