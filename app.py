import streamlit as st
import streamlit.components.v1
import pandas as pd
import numpy as np
import time
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plotly.express as px
import plotly.graph_objects as go

from utils.blast_utils import perform_blast_search, extract_pdb_ids_from_blast, create_multiple_alignment
from utils.pdb_utils import download_pdb_structure, get_structure_info
from utils.visualization import create_structure_viewer, create_alignment_plot, create_msa_visualization, create_conservation_plot
from utils.structure_comparison import (create_structure_comparison_viewer, calculate_structure_alignment, 
                                       create_comparison_metrics_table, create_rmsd_heatmap, 
                                       perform_pairwise_alignment)

# Page configuration
st.set_page_config(
    page_title="BLAST Protein Search & 3D Visualization",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'blast_results' not in st.session_state:
    st.session_state.blast_results = None
if 'pdb_hits' not in st.session_state:
    st.session_state.pdb_hits = []
if 'downloaded_structures' not in st.session_state:
    st.session_state.downloaded_structures = []
if 'sequence_record' not in st.session_state:
    st.session_state.sequence_record = None

# Main title
st.title("🧬 BLAST Protein Search & 3D Visualization")
st.markdown("Search for similar protein structures in the PDB database and visualize them in 3D")

# Sidebar for parameters
st.sidebar.title("Search Parameters")
expect_threshold = st.sidebar.number_input("E-value threshold", value=0.001, format="%.6f", min_value=0.000001)
hitlist_size = st.sidebar.slider("Maximum hits", min_value=5, max_value=50, value=20)
max_downloads = st.sidebar.slider("Max structures to download", min_value=1, max_value=10, value=5)

# Main content tabs
tab1, tab2, tab3, tab4 = st.tabs(["📝 Sequence Input", "🔍 BLAST Results", "🧪 3D Structures", "📊 Analysis & Export"])

with tab1:
    st.header("Protein Sequence Input")
    
    # Input method selection
    input_method = st.radio("Choose input method:", ["Paste sequence", "Upload FASTA file"])
    
    sequence_record = None
    
    if input_method == "Paste sequence":
        col1, col2 = st.columns([2, 1])
        
        with col1:
            sequence_text = st.text_area(
                "Enter protein sequence (amino acids):",
                height=150,
                help="Paste your protein sequence here. Only amino acid sequences are supported."
            )
            
            sequence_id = st.text_input("Sequence ID (optional):", value="query_protein")
            sequence_description = st.text_input("Description (optional):", value="User input protein")
        
        with col2:
            st.info("💡 **Tips:**\n- Enter single-letter amino acid codes\n- Remove any numbers or special characters\n- Minimum length: 10 amino acids")
        
        if sequence_text:
            # Clean sequence
            cleaned_seq = ''.join([c.upper() for c in sequence_text if c.upper() in 'ACDEFGHIKLMNPQRSTVWY'])
            
            if len(cleaned_seq) >= 10:
                sequence_record = SeqRecord(
                    Seq(cleaned_seq), 
                    id=sequence_id, 
                    description=sequence_description
                )
                st.success(f"✅ Sequence loaded: {len(cleaned_seq)} amino acids")
                st.session_state.sequence_record = sequence_record
            else:
                st.error("❌ Sequence too short. Please enter at least 10 amino acids.")
    
    else:  # Upload FASTA file
        uploaded_file = st.file_uploader("Choose a FASTA file", type=['fasta', 'fa', 'fas'])
        
        if uploaded_file is not None:
            try:
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                sequence_record = SeqIO.read(stringio, "fasta")
                st.success(f"✅ File uploaded: {len(sequence_record.seq)} amino acids")
                st.session_state.sequence_record = sequence_record
            except Exception as e:
                st.error(f"❌ Error reading FASTA file: {str(e)}")
    
    # Display sequence info
    if st.session_state.sequence_record:
        with st.expander("Sequence Information", expanded=True):
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Length", len(st.session_state.sequence_record.seq) if st.session_state.sequence_record and st.session_state.sequence_record.seq else 0)
            
            with col2:
                st.metric("ID", st.session_state.sequence_record.id)
            
            with col3:
                st.metric("Description", st.session_state.sequence_record.description[:30] + "..." if len(st.session_state.sequence_record.description) > 30 else st.session_state.sequence_record.description)
            
            seq_str = str(st.session_state.sequence_record.seq) if st.session_state.sequence_record and st.session_state.sequence_record.seq else ""
            st.text_area("Sequence preview:", 
                        value=seq_str[:200] + ("..." if len(seq_str) > 200 else ""),
                        height=100, 
                        disabled=True)

with tab2:
    st.header("BLAST Search Results")
    
    if st.session_state.sequence_record:
        col1, col2 = st.columns([1, 3])
        
        with col1:
            if st.button("🚀 Run BLAST Search", type="primary", use_container_width=True):
                with st.spinner("Performing BLAST search... This may take a few minutes."):
                    try:
                        # Perform BLAST search
                        blast_results = perform_blast_search(
                            st.session_state.sequence_record.seq,
                            expect=expect_threshold,
                            hitlist_size=hitlist_size
                        )
                        
                        if blast_results:
                            st.session_state.blast_results = blast_results
                            
                            # Extract PDB IDs
                            pdb_hits = extract_pdb_ids_from_blast(blast_results, top_n=hitlist_size)
                            st.session_state.pdb_hits = pdb_hits
                            
                            st.success(f"✅ BLAST search completed! Found {len(pdb_hits)} hits")
                        else:
                            st.error("❌ BLAST search failed. Please try again.")
                            
                    except Exception as e:
                        st.error(f"❌ Error during BLAST search: {str(e)}")
        
        with col2:
            if st.session_state.blast_results:
                st.info(f"📊 Query length: {st.session_state.blast_results.seq_len} amino acids")
    
    # Display results
    if st.session_state.pdb_hits:
        st.subheader("BLAST Hits")
        
        # Create results DataFrame
        df_results = pd.DataFrame(st.session_state.pdb_hits)
        
        # Advanced filtering and sorting options
        with st.expander("🔧 Advanced Filtering & Sorting", expanded=False):
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.subheader("Filters")
                # E-value filter
                evalue_filter = st.checkbox("Filter by E-value", value=False)
                if evalue_filter:
                    min_evalue = st.number_input("Min E-value", value=0.0, format="%.2e", min_value=0.0)
                    max_evalue = st.number_input("Max E-value", value=float(df_results['evalue'].max()), format="%.2e", min_value=0.0)
                else:
                    min_evalue, max_evalue = 0.0, float('inf')
                
                # Bit score filter
                bitscore_filter = st.checkbox("Filter by Bit Score", value=False)
                if bitscore_filter:
                    min_bitscore = st.number_input("Min Bit Score", value=float(df_results['bitscore'].min()), min_value=0.0)
                    max_bitscore = st.number_input("Max Bit Score", value=float(df_results['bitscore'].max()), min_value=0.0)
                else:
                    min_bitscore, max_bitscore = 0.0, float('inf')
            
            with col2:
                # Identity filter  
                identity_filter = st.checkbox("Filter by Identity", value=False)
                if identity_filter:
                    min_identity = st.number_input("Min Identity Count", value=int(df_results['identity'].min()), min_value=0)
                    max_identity = st.number_input("Max Identity Count", value=int(df_results['identity'].max()), min_value=0)
                else:
                    min_identity, max_identity = 0, float('inf')
                
                # Alignment length filter
                alignment_filter = st.checkbox("Filter by Alignment Length", value=False)
                if alignment_filter:
                    min_alignment = st.number_input("Min Alignment Length", value=int(df_results['alignment_length'].min()), min_value=1)
                    max_alignment = st.number_input("Max Alignment Length", value=int(df_results['alignment_length'].max()), min_value=1)
                else:
                    min_alignment, max_alignment = 0, float('inf')
            
            with col3:
                st.subheader("Sorting")
                # Sort options
                sort_column = st.selectbox(
                    "Sort by",
                    options=['rank', 'evalue', 'bitscore', 'identity', 'alignment_length', 'pdb_id'],
                    index=1  # Default to E-value
                )
                sort_ascending = st.radio("Sort order", ["Ascending", "Descending"], index=0) == "Ascending"
                
                # PDB ID pattern filter
                pdb_pattern_filter = st.checkbox("Filter by PDB ID pattern", value=False)
                if pdb_pattern_filter:
                    pdb_pattern = st.text_input("PDB ID contains (case-insensitive):", value="")
                else:
                    pdb_pattern = ""
                
                # Apply filters button
                col_apply, col_reset = st.columns(2)
                with col_apply:
                    apply_filters = st.button("🔄 Apply Filters & Sort", type="primary", use_container_width=True)
                with col_reset:
                    reset_filters = st.button("🔄 Reset All", use_container_width=True)
                
                # Reset filters functionality
                if reset_filters:
                    st.session_state.apply_filters = False
                    st.rerun()
        
        # Apply filtering and sorting
        if 'apply_filters' not in st.session_state:
            st.session_state.apply_filters = False
        
        if apply_filters or st.session_state.apply_filters:
            st.session_state.apply_filters = True
            
            # Apply filters
            filtered_df = df_results.copy()
            
            # E-value filter
            if evalue_filter:
                filtered_df = filtered_df[
                    (filtered_df['evalue'] >= min_evalue) & 
                    (filtered_df['evalue'] <= max_evalue)
                ]
            
            # Bit score filter
            if bitscore_filter:
                filtered_df = filtered_df[
                    (filtered_df['bitscore'] >= min_bitscore) & 
                    (filtered_df['bitscore'] <= max_bitscore)
                ]
            
            # Identity filter
            if identity_filter:
                filtered_df = filtered_df[
                    (filtered_df['identity'] >= min_identity) & 
                    (filtered_df['identity'] <= max_identity)
                ]
            
            # Alignment length filter
            if alignment_filter:
                filtered_df = filtered_df[
                    (filtered_df['alignment_length'] >= min_alignment) & 
                    (filtered_df['alignment_length'] <= max_alignment)
                ]
            
            # PDB ID pattern filter
            if pdb_pattern_filter and pdb_pattern:
                # Ensure we have a DataFrame and convert PDB column to string if needed
                if len(filtered_df) > 0:
                    filtered_df = filtered_df[
                        filtered_df['pdb_id'].astype(str).str.contains(pdb_pattern, case=False, na=False)
                    ]
            
            # Apply sorting
            if len(filtered_df) > 0 and sort_column in filtered_df.columns:
                filtered_df = filtered_df.sort_values(by=sort_column, ascending=sort_ascending).reset_index(drop=True)
            
            # Update rank based on new order
            filtered_df['filtered_rank'] = range(1, len(filtered_df) + 1)
            
            df_results = filtered_df
            
            # Show filter results
            if len(filtered_df) < len(st.session_state.pdb_hits):
                st.info(f"📊 Showing {len(filtered_df)} of {len(st.session_state.pdb_hits)} hits after filtering")
            else:
                st.info(f"📊 Showing all {len(filtered_df)} hits (sorted by {sort_column})")
        
        # Display interactive table  
        # Ensure df_results is defined, fallback to original data if not
        if 'df_results' not in locals() or df_results is None:
            df_results = pd.DataFrame(st.session_state.pdb_hits)
            
        display_columns = ['rank', 'pdb_id', 'evalue', 'bitscore', 'identity', 'alignment_length']
        if 'filtered_rank' in df_results.columns:
            display_columns[0] = 'filtered_rank'
        
        st.dataframe(
            df_results[display_columns],
            use_container_width=True,
            hide_index=True,
            column_config={
                "rank": st.column_config.NumberColumn("Rank", width="small"),
                "filtered_rank": st.column_config.NumberColumn("Rank", width="small"),
                "pdb_id": st.column_config.TextColumn("PDB ID", width="small"),
                "evalue": st.column_config.NumberColumn("E-value", format="%.2e"),
                "bitscore": st.column_config.NumberColumn("Bit Score", format="%.1f"),
                "identity": st.column_config.NumberColumn("Identity", width="small"),
                "alignment_length": st.column_config.NumberColumn("Alignment Length", width="small")
            }
        )
        
        # Results visualization
        st.subheader("Results Visualization")
        
        viz_col1, viz_col2 = st.columns(2)
        
        with viz_col1:
            # E-value distribution plot
            st.subheader("E-value Distribution")
            # Prepare data for log scale (clip zeros to avoid log issues)
            df_plot = df_results.copy()
            df_plot['evalue'] = pd.to_numeric(df_plot['evalue'], errors='coerce')
            df_plot['evalue'] = np.where(df_plot['evalue'] <= 0, 1e-300, df_plot['evalue'])
            
            fig = px.bar(df_plot, x='pdb_id', y='evalue', 
                         title=f"E-values of {len(df_plot)} BLAST Hits",
                         labels={'evalue': 'E-value', 'pdb_id': 'PDB ID'})
            fig.update_yaxes(type="log")  # Correct plotly syntax for log scale
            st.plotly_chart(fig, use_container_width=True)
        
        with viz_col2:
            # Identity vs Bit Score scatter
            st.subheader("Identity vs Bit Score")
            fig_scatter = px.scatter(df_results, x='identity', y='bitscore',
                                   hover_data=['pdb_id', 'evalue'],
                                   color='evalue',
                                   size='alignment_length',
                                   title=f"Identity vs Bit Score ({len(df_results)} hits)",
                                   labels={'identity': 'Identity Count', 'bitscore': 'Bit Score', 'evalue': 'E-value'})
            fig_scatter.update_layout(coloraxis_colorbar=dict(title="E-value"))
            st.plotly_chart(fig_scatter, use_container_width=True)
        
        # Multiple Sequence Alignment section
        st.subheader("🧬 Multiple Sequence Alignment")
        
        if st.session_state.sequence_record:
            with st.expander("MSA Visualization Options", expanded=False):
                col1, col2 = st.columns(2)
                
                with col1:
                    max_sequences = st.slider("Max sequences in alignment", min_value=2, max_value=20, value=8)
                    show_conservation = st.checkbox("Show conservation scores", value=True)
                
                with col2:
                    create_msa = st.button("🔬 Create Multiple Alignment", type="secondary", use_container_width=True)
            
            # Create and display MSA
            if create_msa or 'msa_data' in st.session_state:
                if create_msa:
                    with st.spinner("Creating multiple sequence alignment..."):
                        msa_data = create_multiple_alignment(
                            st.session_state.sequence_record.seq,
                            df_results.to_dict('records'),
                            max_sequences=max_sequences
                        )
                        st.session_state.msa_data = msa_data
                else:
                    msa_data = st.session_state.get('msa_data')
                
                if msa_data:
                    st.success(f"✅ Alignment created with {len(msa_data['sequences'])} sequences")
                    
                    # Display MSA visualization
                    col1, col2 = st.columns([3, 1])
                    
                    with col1:
                        st.subheader("Alignment View")
                        fig_msa = create_msa_visualization(msa_data, show_conservation=show_conservation)
                        st.plotly_chart(fig_msa, use_container_width=True)
                    
                    with col2:
                        st.subheader("Statistics")
                        st.metric("Sequences", len(msa_data['sequences']))
                        st.metric("Alignment Length", msa_data['alignment_length'])
                        
                        if 'conservation' in msa_data:
                            avg_conservation = np.mean(msa_data['conservation'])
                            st.metric("Avg Conservation", f"{avg_conservation:.3f}")
                            
                            high_conservation = sum(1 for c in msa_data['conservation'] if c >= 0.8)
                            st.metric("Highly Conserved Positions", high_conservation)
                    
                    # Conservation profile
                    if show_conservation and 'conservation' in msa_data:
                        st.subheader("Conservation Profile")
                        fig_conservation = create_conservation_plot(msa_data)
                        st.plotly_chart(fig_conservation, use_container_width=True)
                    
                    # MSA Export options
                    with st.expander("Export MSA Data", expanded=False):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            # Export alignment as FASTA
                            fasta_content = ""
                            for seq in msa_data['sequences']:
                                fasta_content += f">{seq['id']} {seq['description']}\n{seq['sequence']}\n"
                            
                            st.download_button(
                                label="📄 Download FASTA",
                                data=fasta_content,
                                file_name=f"msa_alignment_{time.strftime('%Y%m%d_%H%M%S')}.fasta",
                                mime="text/plain",
                                use_container_width=True
                            )
                        
                        with col2:
                            # Export conservation scores
                            if 'conservation' in msa_data:
                                conservation_csv = "Position,Conservation_Score\n"
                                for i, score in enumerate(msa_data['conservation'], 1):
                                    conservation_csv += f"{i},{score:.4f}\n"
                                
                                st.download_button(
                                    label="📊 Download Conservation",
                                    data=conservation_csv,
                                    file_name=f"msa_conservation_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                    mime="text/csv",
                                    use_container_width=True
                                )
                else:
                    st.error("❌ Failed to create multiple sequence alignment")
        else:
            st.info("👆 Load a query sequence first to create multiple sequence alignment")
        
    else:
        st.info("👆 Run a BLAST search to see results here")

with tab3:
    st.header("3D Structure Visualization")
    
    if st.session_state.pdb_hits:
        col1, col2 = st.columns([1, 3])
        
        with col1:
            if st.button("📥 Download Structures", type="primary", use_container_width=True):
                with st.spinner("Downloading PDB structures..."):
                    downloaded = []
                    
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    for i, pdb_info in enumerate(st.session_state.pdb_hits[:max_downloads]):
                        pdb_id = pdb_info['pdb_id']
                        status_text.text(f"Downloading {pdb_id}...")
                        
                        file_path = download_pdb_structure(pdb_id)
                        if file_path:
                            pdb_info['file_path'] = file_path
                            downloaded.append(pdb_info)
                        
                        progress_bar.progress((i + 1) / min(max_downloads, len(st.session_state.pdb_hits)))
                        time.sleep(0.5)
                    
                    st.session_state.downloaded_structures = downloaded
                    status_text.text(f"✅ Downloaded {len(downloaded)} structures")
                    st.success(f"Successfully downloaded {len(downloaded)} structures!")
        
        # Structure selection and visualization
        if st.session_state.downloaded_structures:
            st.subheader("Structure Viewer")
            
            # Structure selector
            structure_options = {f"{s['pdb_id']} (E-value: {s['evalue']:.2e})": s for s in st.session_state.downloaded_structures}
            selected_structure_key = st.selectbox("Select structure to visualize:", list(structure_options.keys()))
            
            if selected_structure_key:
                selected_structure = structure_options[selected_structure_key]
                
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    # 3D visualization
                    st.subheader(f"3D Structure: {selected_structure['pdb_id']}")
                    
                    try:
                        viewer = create_structure_viewer(selected_structure['file_path'])
                        streamlit.components.v1.html(viewer, height=500)
                    except Exception as e:
                        st.error(f"Error creating 3D viewer: {str(e)}")
                
                with col2:
                    # Structure information
                    st.subheader("Structure Info")
                    info_data = {
                        "PDB ID": selected_structure['pdb_id'],
                        "E-value": f"{selected_structure['evalue']:.2e}",
                        "Bit Score": f"{selected_structure['bitscore']:.1f}",
                        "Identity": selected_structure['identity'],
                        "Alignment Length": selected_structure['alignment_length'],
                        "Query Range": f"{selected_structure['query_start']}-{selected_structure['query_end']}",
                        "Hit Range": f"{selected_structure['hit_start']}-{selected_structure['hit_end']}"
                    }
                    
                    for key, value in info_data.items():
                        st.metric(key, value)
                    
                    # Download PDB file
                    if os.path.exists(selected_structure['file_path']):
                        with open(selected_structure['file_path'], 'r') as f:
                            pdb_content = f.read()
                        
                        st.download_button(
                            label="📥 Download PDB File",
                            data=pdb_content,
                            file_name=f"{selected_structure['pdb_id']}.pdb",
                            mime="text/plain",
                            use_container_width=True
                        )
        
        else:
            st.info("👆 Download structures to visualize them in 3D")
        
        # Structure Comparison Section
        if st.session_state.downloaded_structures and len(st.session_state.downloaded_structures) >= 2:
            st.divider()
            st.subheader("🔬 Structure Comparison")
            
            comparison_tab1, comparison_tab2 = st.tabs(["Overlay Viewer", "Alignment Metrics"])
            
            with comparison_tab1:
                st.write("Compare multiple structures by overlaying them in a single 3D viewer:")
                
                # Structure selection for comparison
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write("**Select structures to compare:**")
                    structure_selection = {}
                    colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan']
                    
                    for i, structure in enumerate(st.session_state.downloaded_structures[:6]):  # Limit to 6 structures
                        color = colors[i % len(colors)]
                        selected = st.checkbox(
                            f"{structure['pdb_id']} (E-val: {structure['evalue']:.2e})",
                            value=i < 3,  # Select first 3 by default
                            key=f"compare_{structure['pdb_id']}"
                        )
                        if selected:
                            structure_selection[structure['pdb_id']] = {
                                'structure': structure,
                                'color': color
                            }
                    
                    if len(structure_selection) >= 2:
                        if st.button("🔬 Generate Comparison", type="primary", use_container_width=True):
                            st.session_state.comparison_data = structure_selection
                
                with col2:
                    # Display comparison viewer
                    if hasattr(st.session_state, 'comparison_data') and st.session_state.comparison_data:
                        try:
                            # Prepare data for comparison viewer
                            structure_data_list = []
                            for pdb_id, data in st.session_state.comparison_data.items():
                                structure_data_list.append({
                                    'pdb_path': data['structure']['file_path'],
                                    'pdb_id': pdb_id,
                                    'color': data['color'],
                                    'label': f"{pdb_id}"
                                })
                            
                            # Create comparison viewer
                            comparison_viewer = create_structure_comparison_viewer(structure_data_list)
                            st.subheader(f"Comparing {len(structure_data_list)} Structures")
                            streamlit.components.v1.html(comparison_viewer, height=600)
                            
                            # Legend
                            st.write("**Structure Colors:**")
                            legend_cols = st.columns(len(structure_data_list))
                            for i, struct_data in enumerate(structure_data_list):
                                with legend_cols[i]:
                                    st.markdown(f"🔵 **{struct_data['pdb_id']}** - {struct_data['color']}")
                                    
                        except Exception as e:
                            st.error(f"Error creating comparison viewer: {str(e)}")
                    else:
                        st.info("👈 Select at least 2 structures to compare")
            
            with comparison_tab2:
                st.write("Structural alignment metrics and similarity analysis:")
                
                if len(st.session_state.downloaded_structures) >= 2:
                    col1, col2 = st.columns([1, 1])
                    
                    with col1:
                        chain_id = st.selectbox("Select chain for alignment:", ['A', 'B', 'C', 'D'], index=0)
                        
                        if st.button("⚖️ Calculate Pairwise Alignments", type="primary", use_container_width=True):
                            with st.spinner("Calculating structural alignments..."):
                                # Prepare structure paths and labels
                                structure_paths = [s['file_path'] for s in st.session_state.downloaded_structures]
                                structure_labels = [s['pdb_id'] for s in st.session_state.downloaded_structures]
                                
                                # Perform pairwise alignments
                                rmsd_matrix, alignment_results = perform_pairwise_alignment(
                                    structure_paths, structure_labels, chain_id
                                )
                                
                                # Store results in session state
                                st.session_state.alignment_results = alignment_results
                                st.session_state.rmsd_matrix = rmsd_matrix
                                st.session_state.structure_labels = structure_labels
                    
                    with col2:
                        # Display results if available
                        if hasattr(st.session_state, 'alignment_results') and st.session_state.alignment_results:
                            # Metrics table
                            metrics_df = create_comparison_metrics_table(st.session_state.alignment_results)
                            if not metrics_df.empty:
                                st.subheader("Alignment Metrics")
                                st.dataframe(metrics_df, use_container_width=True)
                                
                                # Download metrics
                                csv_data = metrics_df.to_csv(index=False)
                                st.download_button(
                                    label="📊 Download Metrics",
                                    data=csv_data,
                                    file_name=f"structure_comparison_metrics_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                                    mime="text/csv",
                                    use_container_width=True
                                )
                    
                    # RMSD Heatmap
                    if hasattr(st.session_state, 'rmsd_matrix') and st.session_state.rmsd_matrix.size > 0:
                        st.subheader("RMSD Similarity Matrix")
                        heatmap_fig = create_rmsd_heatmap(st.session_state.rmsd_matrix, st.session_state.structure_labels)
                        st.plotly_chart(heatmap_fig, use_container_width=True)
                else:
                    st.info("Need at least 2 downloaded structures for comparison analysis")
    
    else:
        st.info("🔍 Perform a BLAST search first to get structures for visualization")

with tab4:
    st.header("Analysis & Export")
    
    if st.session_state.pdb_hits:
        # Summary statistics
        st.subheader("📊 Summary Statistics")
        
        # Use filtered results if available, otherwise use all results
        if st.session_state.pdb_hits:
            if 'df_results' in locals() and df_results is not None and len(df_results) < len(st.session_state.pdb_hits):
                analysis_df = df_results  # Use filtered results
                st.info(f"📊 Statistics based on {len(analysis_df)} filtered results")
            else:
                analysis_df = pd.DataFrame(st.session_state.pdb_hits)  # Use all results
        else:
            analysis_df = pd.DataFrame()
        
        if not analysis_df.empty:
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Hits", len(analysis_df))
            
            with col2:
                st.metric("Avg E-value", f"{analysis_df['evalue'].mean():.2e}")
            
            with col3:
                st.metric("Avg Identity", f"{analysis_df['identity'].mean():.1f}")
            
            with col4:
                st.metric("Avg Bit Score", f"{analysis_df['bitscore'].mean():.1f}")
        else:
            st.warning("No data available for statistics.")
        
        # Detailed statistics
        if not analysis_df.empty:
            st.subheader("Detailed Statistics")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Identity distribution
                fig_identity = px.histogram(analysis_df, x='identity', 
                                          title="Identity Distribution",
                                          labels={'identity': 'Identity', 'count': 'Frequency'})
                st.plotly_chart(fig_identity, use_container_width=True)
            
            with col2:
                # Bit score vs E-value
                df_scatter = analysis_df.copy()
                df_scatter['evalue'] = pd.to_numeric(df_scatter['evalue'], errors='coerce')
                df_scatter['evalue'] = np.where(df_scatter['evalue'] <= 0, 1e-300, df_scatter['evalue'])
                
                fig_scatter = px.scatter(df_scatter, x='evalue', y='bitscore',
                                       hover_data=['pdb_id', 'identity'],
                                       title="Bit Score vs E-value",
                                       labels={'evalue': 'E-value', 'bitscore': 'Bit Score'})
                fig_scatter.update_xaxes(type="log")
                st.plotly_chart(fig_scatter, use_container_width=True)
        
            # Alignment visualization
            st.subheader("Alignment Region Analysis")
            
            if st.session_state.sequence_record and st.session_state.sequence_record.seq:
                fig_alignment = create_alignment_plot(analysis_df, len(st.session_state.sequence_record.seq))
                st.plotly_chart(fig_alignment, use_container_width=True)
        
            # Export options
            st.subheader("📤 Export Results")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # Export CSV
                csv_data = analysis_df.to_csv(index=False)
                st.download_button(
                    label="📋 Download CSV",
                    data=csv_data,
                    file_name=f"blast_results_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col2:
                # Export JSON
                json_data = analysis_df.to_json(orient='records', indent=2)
                if json_data:
                    st.download_button(
                        label="📄 Download JSON",
                        data=json_data.encode('utf-8'),
                        file_name=f"blast_results_{time.strftime('%Y%m%d_%H%M%S')}.json",
                        mime="application/json",
                        use_container_width=True
                    )
                else:
                    st.error("No data available for JSON export")
        
            with col3:
                # Export summary report
                if st.session_state.sequence_record and not analysis_df.empty:
                    report = f"""BLAST Search Report
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Query Information:
- ID: {st.session_state.sequence_record.id}
- Description: {st.session_state.sequence_record.description}
- Length: {len(st.session_state.sequence_record.seq) if st.session_state.sequence_record.seq else 0} amino acids

Search Parameters:
- E-value threshold: {expect_threshold}
- Maximum hits: {hitlist_size}

Results Summary:
- Total hits found: {len(analysis_df)}
- Average E-value: {analysis_df['evalue'].mean():.2e}
- Average identity: {analysis_df['identity'].mean():.1f}
- Average bit score: {analysis_df['bitscore'].mean():.1f}

Top 5 Hits:
"""
                    for i, (_, hit) in enumerate(analysis_df.head().iterrows(), 1):
                        report += f"{i}. {hit['pdb_id']} - E-value: {hit['evalue']:.2e}, Identity: {hit['identity']}\n"
                    
                    st.download_button(
                        label="📝 Download Report",
                        data=report,
                        file_name=f"blast_report_{time.strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain",
                        use_container_width=True
                    )
        else:
            st.warning("No results available for analysis and export.")
    
    else:
        st.info("🔍 Perform a BLAST search to see analysis and export options")

# Footer
st.markdown("---")
st.markdown("Built with Streamlit • Powered by Biopython & PDB")
