import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="Differential Expression Analysis")
st.title("Differential Expression Analysis: Group Comparison")

# Check if group datasets exist in session state
if 'group_datasets' not in st.session_state or not st.session_state.group_datasets:
    st.warning("Please load and prepare datasets on the Load Data page first!")
    st.stop()

# Get available groups
groups = st.session_state.group_datasets['groups']
datasets = st.session_state.group_datasets['datasets']

if len(groups) < 2:
    st.error("Need at least 2 groups for comparison")
    st.stop()

# Create two columns for group selection
col1, col2 = st.columns(2)

with col1:
    group1 = st.selectbox(
        "Select Group 1 (Control):",
        groups,
        index=0,
        key="group1_selector"
    )

with col2:
    # Filter out group1 from available options
    available_groups = [g for g in groups if g != group1]
    group2 = st.selectbox(
        "Select Group 2 (Case):",
        available_groups,
        index=0 if len(available_groups) > 0 else None,
        key="group2_selector"
    )

# Get data for selected groups
df_group1 = datasets[group1]
df_group2 = datasets[group2]

# Calculate differential expression statistics


def calculate_de_stats(control, case):
    # Calculate mean expression for each group
    mean_control = control.mean(axis=1)
    mean_case = case.mean(axis=1)

    # Calculate fold change (log2)
    fold_change = mean_case - mean_control

    # Calculate p-values (using t-test)
    from scipy.stats import ttest_ind
    p_values = []
    for idx in control.index:
        t, p = ttest_ind(case.loc[idx], control.loc[idx], equal_var=False)
        p_values.append(p)

    # Create results dataframe
    results = pd.DataFrame({
        'gene': control.index,
        'log2_fold_change': fold_change,
        'p_value': p_values,
        '-log10_pvalue': -np.log10(p_values)
    })

    return results


# Добавляем слайдер для выбора количества отображаемых генов
top_genes = st.slider(
    "Number of top significant genes to label:",
    min_value=0,
    max_value=50,
    value=10,
    key="top_genes_slider"
)

if st.button("Run Differential Expression Analysis"):
    with st.spinner("Calculating differential expression..."):
        de_results = calculate_de_stats(df_group1, df_group2)
    
        
        # Сортируем гены по значимости (по p-value)
        de_results_sorted = de_results.sort_values('p_value')
        
        # Создаем колонку для отображения только выбранного количества генов
        de_results_sorted['show_text'] = False
        de_results_sorted.iloc[:top_genes, de_results_sorted.columns.get_loc('show_text')] = True

        # Create volcano plot
        fig = px.scatter(
            de_results_sorted,
            x='log2_fold_change',
            y='-log10_pvalue',
            text='gene',
            hover_name='gene',
            labels={
                'log2_fold_change': 'log2(Fold Change)',
                '-log10_pvalue': '-log10(p-value)'
            },
            title=f"Volcano Plot: {group1} vs {group2}"
        )
        
        # Обновляем отображение текста только для выбранных генов
        fig.update_traces(
            textposition='top center',
            textfont_size=10,
            texttemplate='%{text}',
            selector=dict(mode='markers+text')
        )
        
        # Показываем текст только для генов с show_text=True
        fig.for_each_trace(
            lambda trace: trace.update(textfont_color=['white' if show else 'rgba(0,0,0,0)' 
                                      for show in de_results_sorted['show_text']])
        )

        # Add significance thresholds (optional)
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="red")
        fig.add_vline(x=1, line_dash="dash", line_color="blue")
        fig.add_vline(x=-1, line_dash="dash", line_color="blue")

        st.plotly_chart(fig, use_container_width=True)

        # Show results table
        st.subheader("Differential Expression Results")
        st.dataframe(de_results_sorted.sort_values('p_value'))
        
        # st.dataframe(de_results_sorted.sort_values('p_value')[:top_genes])


        # Create the heatmap
        fig = plt.figure()  # Adjust the figure size
        filtered = st.session_state.group_datasets['avg'].loc[:, st.session_state.group_datasets['avg'].columns.isin(de_results_sorted.sort_values('p_value')[:top_genes]['gene'])]
        sns.heatmap(filtered, annot=True, fmt='.1f', linewidths=0.5)

        # Add labels and title
        plt.title('Gene Expression Heatmap')
        plt.xlabel('Gene')
        plt.ylabel('Group')

        # Show the plot
        st.subheader("Heatmap")
        st.pyplot(fig)