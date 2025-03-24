import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_heatmap = st.file_uploader(label='Upload your dataset:', key=0)

if df_heatmap:
    df_heatmap = pd.read_excel(df_heatmap)
    df_pivot = df_heatmap.pivot_table(index='Group', columns='Gene', values='Expression')
    df_pivot.head()

    # Calculate the range (max - min) for each gene
    ranges = df_pivot.max() - df_pivot.min()

    # Filter columns where the range is at least 1.0
    filtered_genes = ranges[ranges >= 6.0].index

    # Keep only the filtered genes in the DataFrame
    df_filtered = df_pivot[filtered_genes]

    print(df_filtered)

    # Create the heatmap
    fig = plt.figure(figsize=(50, 6))  # Adjust the figure size
    sns.heatmap(df_filtered, annot=True, fmt='.1f', linewidths=0.5)

    # Add labels and title
    plt.title('Gene Expression Heatmap')
    plt.xlabel('Gene')
    plt.ylabel('Group')

    # Show the plot
    st.pyplot(fig)