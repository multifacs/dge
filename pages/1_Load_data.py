import streamlit as st
import os
import pandas as pd
import subprocess
from pathlib import Path
from datetime import datetime

# Config
DATA_DIR = "data"
Path(DATA_DIR).mkdir(exist_ok=True)

st.title("GEO Data File Processor")

# Get files in data directory
@st.cache_data(ttl=5)  # Refresh every 5 seconds
def get_files(extension='.gz'):
    files = []
    for item in os.listdir(DATA_DIR):
        item_path = os.path.join(DATA_DIR, item)
        if os.path.isfile(item_path) and item.endswith(extension):
            files.append({
                "Filename": item,
                "Size (KB)": round(os.path.getsize(item_path)/1024, 2),
                "Last Modified": datetime.fromtimestamp(os.path.getmtime(item_path))
            })
    return pd.DataFrame(files).sort_values("Last Modified", ascending=False)

# First file selection (GZ files)
files_df = get_files('.gz')

if files_df.empty:
    st.warning("No GZ files found in data folder")
    st.stop()

# Display selectable table for GZ files
st.subheader("Available GEO datasets")
gz_selection = st.dataframe(
    files_df,
    key="gz_file_selector",
    on_select="rerun",
    selection_mode="single-row",
    use_container_width=True
)

# Process selected GEO file
if gz_selection['selection']['rows']:
    selected_row = gz_selection['selection']['rows'][0]
    selected_file = files_df.iloc[selected_row]["Filename"]
    input_path = os.path.join(DATA_DIR, selected_file)
    
    name = input_path.split("\\")[1].split("_")[0]

    with st.spinner(f"Processing {selected_file}..."):
        try:
            if not os.path.exists(f"data//{name}//{name}_expr.csv"):
                # Run R script
                result = subprocess.run(
                    ["Rscript", "process_geo.R", input_path, name],
                    capture_output=True,
                    text=True,
                    check=True
                )
                st.success(f"Successfully processed {selected_file}!")
            
            if os.path.exists(f"data//{name}//{name}_expr.csv"):
                # Load and process expression data
                st.subheader("Expression matrix")
                gene_names_df = pd.read_excel("data/conv_table.xlsx", index_col=0).rename_axis('id')
                expr_df = pd.read_csv(f"data//{name}//{name}_expr.csv", index_col=0).rename_axis('id')
                expr_df.index = expr_df.index.map(lambda x: gene_names_df.loc[x, 'symbol'])
                expr_df.rename_axis('symbol', inplace=True)
                st.dataframe(expr_df)
                
                # Show phenotype data
                st.subheader("Phenotype data")
                phen_df = pd.read_csv(f"data//{name}//{name}_phen.csv", index_col=0).rename_axis('id')
                st.dataframe(phen_df)
                
                # Second file selection (TXT gene lists)
                st.markdown("---")
                st.subheader("Gene list selection")
                txt_files_df = get_files('.txt')
                
                if txt_files_df.empty:
                    st.warning("No TXT gene list files found in data folder")
                else:
                    # Display selectable table for TXT files
                    txt_selection = st.dataframe(
                        txt_files_df,
                        key="txt_file_selector",
                        on_select="rerun",
                        selection_mode="single-row",
                        use_container_width=True
                    )
                    
                    if txt_selection['selection']['rows']:
                        selected_txt_row = txt_selection['selection']['rows'][0]
                        selected_txt_file = txt_files_df.iloc[selected_txt_row]["Filename"]
                        txt_path = os.path.join(DATA_DIR, selected_txt_file)
                        
                        # Read gene list
                        with open(txt_path, 'r') as f:
                            gene_list = [line.strip() for line in f.readlines() if line.strip()]
                        
                        # Filter expression matrix
                        filtered_expr = expr_df[expr_df.index.isin(gene_list)]
                        
                        st.subheader(f"Filtered Expression Matrix (by {selected_txt_file})")
                        st.write(f"Genes in list: {len(gene_list)}")
                        st.write(f"Genes found in matrix: {len(filtered_expr)}")
                        st.dataframe(filtered_expr)
                        
                        # Add phenotype column selector for additional filtering
                        if not filtered_expr.empty:
                            st.markdown("---")
                            st.subheader("Additional filtering by phenotype")
                            
                            # Get all phenotype columns (excluding index)
                            phen_columns = [col for col in phen_df.columns if col != phen_df.index.name]
                            
                            if phen_columns:
                                selected_phen_column = st.selectbox(
                                    "Select phenotype column to filter samples:",
                                    phen_columns,
                                    key="phen_column_selector"
                                )
                                
                                # Get unique values in selected column
                                phen_values = phen_df[selected_phen_column].unique()
                                selected_values = st.multiselect(
                                    f"Select {selected_phen_column} values to keep:",
                                    phen_values,
                                    default=phen_values,  # Show all by default
                                    key="phen_value_selector"
                                )
                                
                                if selected_values:
                                    # Get sample IDs matching the selected phenotype values
                                    filtered_samples = phen_df[phen_df[selected_phen_column].isin(selected_values)].index.tolist()
                                    
                                    # Filter expression matrix columns (keep gene column + only matching samples)
                                    final_filtered = filtered_expr[filtered_samples]
                                    
                                    st.subheader("Final Filtered Expression Matrix")
                                    st.write(f"Samples matching criteria: {len(filtered_samples)}")
                                    st.dataframe(final_filtered)
                                    
                                # Add button to create datasets by groups
                                if st.button("Создать датасеты по группам", key="create_group_datasets"):
                                    if 'group_datasets' not in st.session_state:
                                        st.session_state.group_datasets = {}
                                    
                                    # Get all unique groups in selected column
                                    unique_groups = phen_df[selected_phen_column].unique()
                                    
                                    # Create datasets for each group
                                    group_datasets = {}
                                    for group in unique_groups:
                                        # Get sample IDs for this group
                                        group_samples = phen_df[phen_df[selected_phen_column] == group].index.tolist()
                                        
                                        # Filter expression matrix for this group
                                        group_data = filtered_expr[group_samples]
                                        
                                        # Store in dictionary
                                        group_datasets[group] = group_data
                                    
                                    # Save to session state
                                    st.session_state.group_datasets = {
                                        'groups': list(group_datasets.keys()),
                                        'datasets': group_datasets
                                    }
                                    
                                    st.success(f"Создано {len(unique_groups)} датасетов по группам!")

                                # Display created datasets if they exist
                                if 'group_datasets' in st.session_state and st.session_state.group_datasets:
                                    st.markdown("---")
                                    st.subheader("Датасеты по группам")
                                    
                                    tabs = st.tabs([f"Группа {group}" for group in st.session_state.group_datasets['groups']])
                                    
                                    for tab, group in zip(tabs, st.session_state.group_datasets['groups']):
                                        with tab:
                                            dataset = st.session_state.group_datasets['datasets'][group]
                                            st.write(f"Образцов в группе: {len(dataset.columns)}")
                                            st.dataframe(dataset)
                            else:
                                st.warning("No phenotype columns available for filtering")
            
        except subprocess.CalledProcessError as e:
            st.error(f"Processing failed: {e.stderr}")
else:
    st.info("👆 Select a GEO dataset file from the table above to process")

# Add some debug info (can be removed)
with st.expander("Debug Info"):
    st.write("GZ Selection object:", gz_selection)
    st.write("TXT Files DataFrame:", get_files('.txt'))