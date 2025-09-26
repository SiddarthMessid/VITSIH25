"""
3D visualization utilities for protein structures
"""

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import py3Dmol
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.vectors import Vector


def create_structure_viewer(pdb_file_path, width=600, height=400):
    """
    Create 3D structure viewer using py3Dmol
    
    Args:
        pdb_file_path: Path to PDB file
        width: Viewer width (default: 600)
        height: Viewer height (default: 400)
    
    Returns:
        HTML string for embedding the viewer
    """
    try:
        if not os.path.exists(pdb_file_path):
            return "<p>PDB file not found</p>"
        
        # Read PDB file content
        with open(pdb_file_path, 'r') as f:
            pdb_content = f.read()
        
        # Create py3Dmol viewer
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(pdb_content, 'pdb')
        
        # Set cartoon style for proteins
        viewer.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Add surface representation
        viewer.addSurface(py3Dmol.VDW, {'opacity': 0.3, 'color': 'white'})
        
        # Center and zoom
        viewer.zoomTo()
        
        # Return HTML
        return viewer._make_html()
        
    except Exception as e:
        return f"<p>Error creating 3D viewer: {str(e)}</p>"


def create_alignment_plot(df_results, query_length):
    """
    Create alignment region visualization showing where each BLAST hit aligns to the query sequence
    
    Args:
        df_results: DataFrame with BLAST results containing query_start, query_end, pdb_id, evalue
        query_length: Length of query sequence
    
    Returns:
        Plotly figure object
    """
    try:
        fig = go.Figure()
        
        # Add query sequence as reference bar at the bottom
        fig.add_trace(go.Bar(
            x=[query_length],
            y=["Query"],
            orientation='h',
            name="Query Sequence",
            marker=dict(color="lightgray"),
            text=f"Query ({query_length} AA)",
            textposition="middle center"
        ))
        
        # Add alignment regions for each hit
        colors = px.colors.qualitative.Set1
        y_positions = []
        hover_text = []
        x_starts = []
        x_lengths = []
        bar_colors = []
        
        for i, (_, hit) in enumerate(df_results.iterrows()):
            color = colors[i % len(colors)]
            y_pos = f"{hit['pdb_id']}"
            
            # Calculate alignment length and position
            alignment_length = hit['query_end'] - hit['query_start'] + 1
            
            y_positions.append(y_pos)
            x_starts.append(hit['query_start'])
            x_lengths.append(alignment_length)
            bar_colors.append(color)
            
            # Create hover text with detailed info
            hover_info = (f"PDB: {hit['pdb_id']}<br>"
                         f"Query: {hit['query_start']}-{hit['query_end']}<br>"
                         f"Identity: {hit['identity']}<br>"
                         f"E-value: {hit['evalue']:.2e}<br>"
                         f"Bit Score: {hit['bitscore']:.1f}")
            hover_text.append(hover_info)
        
        # Add alignment bars
        if y_positions:
            fig.add_trace(go.Bar(
                x=x_lengths,
                y=y_positions,
                base=x_starts,
                orientation='h',
                name="Alignments",
                marker=dict(color=bar_colors),
                hovertemplate="%{hovertext}<extra></extra>",
                hovertext=hover_text
            ))
        
        # Update layout
        fig.update_layout(
            title="Alignment Regions on Query Sequence",
            xaxis_title="Amino Acid Position",
            yaxis_title="BLAST Hits",
            xaxis=dict(range=[0, query_length * 1.1]),
            height=max(400, (len(df_results) + 1) * 40 + 100),
            showlegend=False,
            plot_bgcolor='white'
        )
        
        # Add reference lines at regular intervals
        for pos in range(0, query_length, max(1, query_length // 10)):
            fig.add_vline(x=pos, line_dash="dash", line_color="lightgray", opacity=0.5)
        
        return fig
        
    except Exception as e:
        # Return empty figure on error
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating alignment plot: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig


def create_msa_visualization(alignment_data, show_conservation=True):
    """
    Create interactive multiple sequence alignment visualization
    
    Args:
        alignment_data: Alignment data from create_multiple_alignment
        show_conservation: Whether to show conservation scores
    
    Returns:
        Plotly figure object
    """
    try:
        if not alignment_data or not alignment_data['sequences']:
            fig = go.Figure()
            fig.add_annotation(
                text="No alignment data available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig
        
        sequences = alignment_data['sequences']
        alignment_length = alignment_data['alignment_length']
        
        # Create data for text-based alignment display
        seq_texts = []
        seq_labels = []
        
        for i, seq in enumerate(sequences):
            sequence = seq['sequence']
            seq_id = seq['id']
            
            # Pad sequence to alignment length
            padded_seq = sequence.ljust(alignment_length, '-')
            seq_texts.append(list(padded_seq))
            
            # Create sequence label with additional info
            if seq_id == 'Query':
                seq_labels.append(f"Query ({len(sequence)} AA)")
            else:
                evalue = seq.get('evalue', 0)
                identity = seq.get('identity', 0)
                seq_labels.append(f"{seq_id} (E: {evalue:.2e}, I: {identity})")
        
        # Create figure with text annotations
        fig = go.Figure()
        
        # Add sequence text
        for i, (seq_text, label) in enumerate(zip(seq_texts, seq_labels)):
            y_pos = len(seq_texts) - i - 1  # Reverse order
            
            # Add sequence as text
            seq_string = ''.join(seq_text)
            fig.add_annotation(
                x=0, y=y_pos,
                text=f"{label}: {seq_string}",
                showarrow=False,
                xanchor="left",
                font=dict(family="monospace", size=10)
            )
        
        # Update layout for text display
        fig.update_layout(
            title="Multiple Sequence Alignment",
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            height=max(300, len(sequences) * 30 + 100),
            plot_bgcolor='white',
            margin=dict(l=20, r=20, t=50, b=20)
        )
        
        return fig
        
    except Exception as e:
        # Return empty figure on error
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating MSA visualization: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig


def create_conservation_plot(alignment_data):
    """
    Create conservation score plot for multiple sequence alignment
    
    Args:
        alignment_data: Alignment data with conservation scores
    
    Returns:
        Plotly figure object
    """
    try:
        if not alignment_data or 'conservation' not in alignment_data:
            fig = go.Figure()
            fig.add_annotation(
                text="No conservation data available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig
        
        conservation = alignment_data['conservation']
        positions = list(range(1, len(conservation) + 1))
        
        fig = go.Figure()
        
        # Add conservation line plot
        fig.add_trace(go.Scatter(
            x=positions,
            y=conservation,
            mode='lines+markers',
            name='Conservation Score',
            line=dict(color='blue', width=2),
            marker=dict(size=4),
            hovertemplate="Position: %{x}<br>Conservation: %{y:.3f}<extra></extra>"
        ))
        
        # Add high conservation threshold line
        fig.add_hline(y=0.8, line_dash="dash", line_color="red", 
                     annotation_text="High Conservation (0.8)")
        
        fig.update_layout(
            title="Sequence Conservation Profile",
            xaxis_title="Alignment Position",
            yaxis_title="Conservation Score",
            yaxis=dict(range=[0, 1]),
            height=300
        )
        
        return fig
        
    except Exception as e:
        # Return empty figure on error
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating conservation plot: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig


def create_identity_heatmap(df_results):
    """
    Create heatmap showing identity scores
    
    Args:
        df_results: DataFrame with BLAST results
    
    Returns:
        Plotly figure object
    """
    try:
        # Prepare data for heatmap
        data = df_results[['pdb_id', 'identity', 'evalue', 'bitscore']].copy()
        data = data.sort_values('evalue')
        
        fig = go.Figure(data=go.Heatmap(
            z=[data['identity'].values],
            x=data['pdb_id'].values,
            y=['Identity'],
            colorscale='RdYlBu_r',
            text=[[f"E: {e:.1e}<br>Score: {s:.1f}" for e, s in zip(data['evalue'], data['bitscore'])]],
            texttemplate="%{text}",
            textfont={"size": 10},
            colorbar=dict(title="Identity Score")
        ))
        
        fig.update_layout(
            title="Identity Scores Across PDB Hits",
            xaxis_title="PDB ID",
            height=200
        )
        
        return fig
        
    except Exception as e:
        # Return empty figure on error
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating heatmap: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig




def create_3d_scatter_plot(df_results):
    """
    Create 3D scatter plot of BLAST results
    
    Args:
        df_results: DataFrame with BLAST results
    
    Returns:
        Plotly figure object
    """
    try:
        fig = px.scatter_3d(
            df_results,
            x='evalue',
            y='bitscore',
            z='identity',
            color='alignment_length',
            hover_data=['pdb_id', 'description'],
            title="3D Analysis of BLAST Results",
            labels={
                'evalue': 'E-value',
                'bitscore': 'Bit Score',
                'identity': 'Identity',
                'alignment_length': 'Alignment Length'
            }
        )
        
        fig.update_layout(
            scene=dict(
                xaxis=dict(type="log", title="E-value (log scale)"),
                yaxis=dict(title="Bit Score"),
                zaxis=dict(title="Identity")
            ),
            height=600
        )
        
        return fig
        
    except Exception as e:
        # Return empty figure on error
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error creating 3D plot: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
