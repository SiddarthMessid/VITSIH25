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
    Create alignment region visualization
    
    Args:
        df_results: DataFrame with BLAST results
        query_length: Length of query sequence
    
    Returns:
        Plotly figure object
    """
    try:
        fig = go.Figure()
        
        # Add query sequence as reference line
        fig.add_shape(
            type="line",
            x0=1, y0=0, x1=query_length, y1=0,
            line=dict(color="black", width=4),
            name="Query Sequence"
        )
        
        # Add alignment regions for each hit
        colors = px.colors.qualitative.Set1
        
        for i, hit in df_results.iterrows():
            color = colors[i % len(colors)]
            
            # Add alignment region
            fig.add_shape(
                type="line",
                x0=hit['query_start'], y0=i+1,
                x1=hit['query_end'], y1=i+1,
                line=dict(color=color, width=6),
            )
            
            # Add text annotation
            fig.add_annotation(
                x=hit['query_end'] + 10,
                y=i+1,
                text=f"{hit['pdb_id']} (E: {hit['evalue']:.1e})",
                showarrow=False,
                xanchor="left",
                font=dict(size=10)
            )
        
        # Update layout
        fig.update_layout(
            title="Alignment Regions on Query Sequence",
            xaxis_title="Amino Acid Position",
            yaxis_title="BLAST Hits",
            xaxis=dict(range=[0, query_length + 100]),
            yaxis=dict(range=[-0.5, len(df_results) + 0.5]),
            height=max(400, len(df_results) * 30 + 100),
            showlegend=False
        )
        
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
