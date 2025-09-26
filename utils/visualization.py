"""
3D visualization utilities for protein structures
"""

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import py3Dmol
import os


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
