import streamlit as st
import pandas as pd
import time
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plotly.express as px
import plotly.graph_objects as go

from utils.blast_utils import perform_blast_search, extract_pdb_ids_from_blast
from utils.pdb_utils import download_pdb_structure, get_structure_info
from utils.visualization import create_structure_viewer, create_alignment_plot

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
                st.metric("Length", len(st.session_state.sequence_record.seq))
            
            with col2:
                st.metric("ID", st.session_state.sequence_record.id)
            
            with col3:
                st.metric("Description", st.session_state.sequence_record.description[:30] + "..." if len(st.session_state.sequence_record.description) > 30 else st.session_state.sequence_record.description)
            
            st.text_area("Sequence preview:", 
                        value=str(st.session_state.sequence_record.seq)[:200] + ("..." if len(st.session_state.sequence_record.seq) > 200 else ""),
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
        st.subheader("Top BLAST Hits")
        
        # Create results DataFrame
        df_results = pd.DataFrame(st.session_state.pdb_hits)
        
        # Display interactive table
        display_columns = ['rank', 'pdb_id', 'evalue', 'bitscore', 'identity', 'alignment_length']
        st.dataframe(
            df_results[display_columns],
            use_container_width=True,
            hide_index=True,
            column_config={
                "rank": st.column_config.NumberColumn("Rank", width="small"),
                "pdb_id": st.column_config.TextColumn("PDB ID", width="small"),
                "evalue": st.column_config.NumberColumn("E-value", format="%.2e"),
                "bitscore": st.column_config.NumberColumn("Bit Score", format="%.1f"),
                "identity": st.column_config.NumberColumn("Identity", width="small"),
                "alignment_length": st.column_config.NumberColumn("Alignment Length", width="small")
            }
        )
        
        # E-value distribution plot
        st.subheader("E-value Distribution")
        fig = px.bar(df_results, x='pdb_id', y='evalue', 
                     title="E-values of BLAST Hits",
                     labels={'evalue': 'E-value', 'pdb_id': 'PDB ID'})
        fig.update_layout(yscale="log")
        st.plotly_chart(fig, use_container_width=True)
        
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
                        st.components.v1.html(viewer, height=500)
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
    
    else:
        st.info("🔍 Perform a BLAST search first to get structures for visualization")

with tab4:
    st.header("Analysis & Export")
    
    if st.session_state.pdb_hits:
        # Summary statistics
        st.subheader("📊 Summary Statistics")
        
        df_results = pd.DataFrame(st.session_state.pdb_hits)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Hits", len(df_results))
        
        with col2:
            st.metric("Avg E-value", f"{df_results['evalue'].mean():.2e}")
        
        with col3:
            st.metric("Avg Identity", f"{df_results['identity'].mean():.1f}")
        
        with col4:
            st.metric("Avg Bit Score", f"{df_results['bitscore'].mean():.1f}")
        
        # Detailed statistics
        st.subheader("Detailed Statistics")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Identity distribution
            fig_identity = px.histogram(df_results, x='identity', 
                                      title="Identity Distribution",
                                      labels={'identity': 'Identity', 'count': 'Frequency'})
            st.plotly_chart(fig_identity, use_container_width=True)
        
        with col2:
            # Bit score vs E-value
            fig_scatter = px.scatter(df_results, x='evalue', y='bitscore',
                                   hover_data=['pdb_id', 'identity'],
                                   title="Bit Score vs E-value",
                                   labels={'evalue': 'E-value', 'bitscore': 'Bit Score'})
            fig_scatter.update_xaxes(type="log")
            st.plotly_chart(fig_scatter, use_container_width=True)
        
        # Alignment visualization
        st.subheader("Alignment Region Analysis")
        
        if st.session_state.sequence_record:
            fig_alignment = create_alignment_plot(df_results, len(st.session_state.sequence_record.seq))
            st.plotly_chart(fig_alignment, use_container_width=True)
        
        # Export options
        st.subheader("📤 Export Results")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # Export CSV
            csv_data = df_results.to_csv(index=False)
            st.download_button(
                label="📋 Download CSV",
                data=csv_data,
                file_name=f"blast_results_{time.strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv",
                use_container_width=True
            )
        
        with col2:
            # Export JSON
            json_data = df_results.to_json(orient='records', indent=2)
            st.download_button(
                label="📄 Download JSON",
                data=json_data,
                file_name=f"blast_results_{time.strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json",
                use_container_width=True
            )
        
        with col3:
            # Export summary report
            if st.session_state.sequence_record:
                report = f"""BLAST Search Report
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Query Information:
- ID: {st.session_state.sequence_record.id}
- Description: {st.session_state.sequence_record.description}
- Length: {len(st.session_state.sequence_record.seq)} amino acids

Search Parameters:
- E-value threshold: {expect_threshold}
- Maximum hits: {hitlist_size}

Results Summary:
- Total hits found: {len(df_results)}
- Average E-value: {df_results['evalue'].mean():.2e}
- Average identity: {df_results['identity'].mean():.1f}
- Average bit score: {df_results['bitscore'].mean():.1f}

Top 5 Hits:
"""
                for i, hit in enumerate(df_results.head().itertuples(), 1):
                    report += f"{i}. {hit.pdb_id} - E-value: {hit.evalue:.2e}, Identity: {hit.identity}\n"
                
                st.download_button(
                    label="📝 Download Report",
                    data=report,
                    file_name=f"blast_report_{time.strftime('%Y%m%d_%H%M%S')}.txt",
                    mime="text/plain",
                    use_container_width=True
                )
    
    else:
        st.info("🔍 Perform a BLAST search to see analysis and export options")

# Footer
st.markdown("---")
st.markdown("Built with Streamlit • Powered by Biopython & PDB")
