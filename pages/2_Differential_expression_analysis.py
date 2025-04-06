import streamlit as st
from streamlit_theme import st_theme
import pandas as pd
import plotly.express as px
from pathlib import Path
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

st.set_page_config(page_title="Дифференциальный анализ экспрессии")
st.title("Дифференциальный анализ экспрессии")


def validate_session_state():
    if 'group_datasets' not in st.session_state or not st.session_state.group_datasets:
        st.warning(
            "Сначала подготовьте датасет на первой странице!")
        st.stop()
    if len(st.session_state.group_datasets['groups']) < 2:
        st.error("Для анализа нужно как минимум 2 группы")
        st.stop()


def select_groups():
    groups = st.session_state.group_datasets['groups']
    col1, col2 = st.columns(2)
    with col1:
        group1 = st.selectbox("Выберите Группу 1 (Контроль):",
                              groups, index=0, key="group1_selector")
    with col2:
        available_groups = [g for g in groups if g != group1]
        group2 = st.selectbox("Выберите Группу 2 (Случай):",
                              available_groups, index=0, key="group2_selector")
    return group1, group2


def calculate_de_stats(control, case):
    mean_control = control.mean(axis=1)
    mean_case = case.mean(axis=1)
    fold_change = mean_case - mean_control
    p_values = [ttest_ind(case.loc[idx], control.loc[idx], equal_var=False)[
        1] for idx in control.index]
    return pd.DataFrame({
        'gene': control.index,
        'log2_fold_change': fold_change,
        'p_value': p_values,
        '-log10_pvalue': -np.log10(p_values)
    }).reset_index(drop=True)


def get_text_color():
    # theme = dict(st_theme())
    # return theme["textColor"]
    # return '#fafafa'
    return st.get_option('theme.textColor')


def plot_volcano(results, group1, group2, top_genes=10, fc_threshold=1, P_VALUE=0.05):
    BLUE = '#36a2eb'
    RED = '#ff6384'
    GREEN = '#4bc0c0'

    results_sorted = results.sort_values('p_value')
    significant = results_sorted[(results_sorted['log2_fold_change'].abs(
    ) >= fc_threshold) & (results_sorted['p_value'] < P_VALUE)]
    top_label = significant.head(top_genes)

    st.session_state.top_genes = top_label['gene'].tolist()

    text_color = get_text_color()

    results_sorted['is_significant'] = (
        (results_sorted['p_value'] < P_VALUE) &
        (results_sorted['log2_fold_change'].abs() >= fc_threshold)
    )

    # Основной слой — все точки, без текста
    fig = px.scatter(
        results_sorted,
        x='log2_fold_change',
        y='-log10_pvalue',
        hover_name='gene',
        labels={
            'log2_fold_change': 'log2(Fold Change)',
            '-log10_pvalue': '-log10(p-value)'
        },
        title=f"Volcano Plot: {group1} vs {group2}",
        color='is_significant',
        color_discrete_map={True: RED, False: BLUE}
    )

    # Добавляем второй слой — только top genes, с подписями
    top_data = results_sorted.loc[top_label.index]
    fig.add_scatter(
        x=top_data['log2_fold_change'],
        y=top_data['-log10_pvalue'],
        mode='text',
        text=top_data['gene'],
        textposition='top center',
        textfont=dict(color=text_color, size=15),
        showlegend=False,
        hoverinfo='skip'
    )

    # Add threshold lines
    fc_pos = fc_threshold
    fc_neg = -fc_threshold
    pval_line = -np.log10(P_VALUE)

    fig.add_hline(y=pval_line, line_dash="dash", line_color=RED)
    fig.add_vline(x=fc_pos, line_dash="dash", line_color=GREEN)
    fig.add_vline(x=fc_neg, line_dash="dash", line_color=GREEN)

    # Add legend annotations
    fig.add_annotation(
        x=fc_pos, y=results_sorted['-log10_pvalue'].max(),
        text="log2(FC) ≥ 1", showarrow=True, arrowhead=1, ax=40, ay=-30,
        font=dict(color=GREEN)
    )
    fig.add_annotation(
        x=fc_neg, y=results_sorted['-log10_pvalue'].max(),
        text="log2(FC) ≤ -1", showarrow=True, arrowhead=1, ax=-40, ay=-30,
        font=dict(color=GREEN)
    )
    fig.add_annotation(
        x=results_sorted['log2_fold_change'].min(), y=pval_line,
        text=f"p-value = {P_VALUE}", showarrow=True, arrowhead=1, ax=0, ay=-40,
        font=dict(color=RED)
    )

    st.plotly_chart(fig, use_container_width=True)


def plot_heatmap(group1, group2, group1_data, group2_data, top_genes, fc_threshold, P_VALUE=0.05, width=4, height=40):
    # Calculate mean expression per group for top genes
    mean_expr = pd.DataFrame({
        group1: group1_data.loc[top_genes].mean(axis=1),
        group2: group2_data.loc[top_genes].mean(axis=1)
    })

    # print(mean_expr.head(10))
    mean_expr.sort_values('id', inplace=True)
    plt.rcParams.update({'font.size': 6})

    fig = plt.figure(figsize=(width, height))
    sns.heatmap(mean_expr, annot=True, fmt='.2f',
                linewidths=0.5, cmap='plasma')
    plt.title(
        f'Средняя экспрессия топ-{len(top_genes)} генов (|FC| ≥ {fc_threshold} & p-value < {P_VALUE})')
    plt.xlabel('Группа')
    plt.ylabel('Ген')
    st.subheader("Тепловая карта Средней экспрессии по группам")
    st.pyplot(fig)

# def plot_heatmap(group1, group2, group1_data, group2_data, top_genes, fc_threshold):
#     mean_expr = pd.DataFrame({
#         group1: group1_data.loc[top_genes].mean(axis=1),
#         group2: group2_data.loc[top_genes].mean(axis=1)
#     })

#     # mean_expr['abs_diff'] = (mean_expr[group1] - mean_expr[group2]).abs()
#     # mean_expr = mean_expr.sort_values(by='abs_diff', ascending=False).drop(columns='abs_diff')
#     mean_expr.sort_values('id', inplace=True)
#     print(mean_expr.head(10))

#     chunk_size = 50
#     num_chunks = int(np.ceil(len(mean_expr) / chunk_size))

#     st.subheader("Тепловая карта средней экспрессии по группам")
#     for i in range(num_chunks):
#         chunk = mean_expr.iloc[i*chunk_size:(i+1)*chunk_size]

#         fig, ax = plt.subplots(figsize=(4, max(2, 0.2 * len(chunk))))
#         sns.heatmap(chunk, annot=True, fmt='.2f', linewidths=0.5, cmap='plasma', ax=ax)
#         ax.set_title(f'Средняя экспрессия генов {i*chunk_size+1}–{min((i+1)*chunk_size, len(mean_expr))}')
#         ax.set_xlabel('Группа')
#         ax.set_ylabel('Ген')
#         st.pyplot(fig)


def main():
    validate_session_state()
    group1, group2 = select_groups()
    df_group1 = st.session_state.group_datasets['datasets'][group1]
    df_group2 = st.session_state.group_datasets['datasets'][group2]
    fc_threshold = 1

    P_VALUE = float(st.text_input("P-value", "0.05"))
    width = float(st.text_input("Heatmap width", "4"))
    height = float(st.text_input("Heatmap height", "40"))

    if st.button("Запустить дифференциальный анализ экспрессии"):
        with st.spinner("Считаем дифференциальную экспрессию..."):
            results = calculate_de_stats(df_group1, df_group2)
            st.session_state.analysis_results = results
            st.subheader("Все результаты дифф. экспрессии")
            st.dataframe(results.sort_values('p_value'))

    if 'analysis_results' in st.session_state:
        results = st.session_state.analysis_results
        sig_genes = results[(results['p_value'] < P_VALUE) & (
            results['log2_fold_change'].abs() >= fc_threshold)]
        max_top_genes = len(sig_genes)

        if max_top_genes > 0:
            top_genes = st.slider("Число помеченных самых значимых генов:", 0, max_top_genes, min(
                10, max_top_genes), key="top_genes_slider")
            plot_volcano(results, group1, group2,
                         top_genes, fc_threshold, P_VALUE)

            if 'top_genes' in st.session_state and st.session_state.top_genes:
                plot_heatmap(group1, group2, df_group1, df_group2,
                             st.session_state.top_genes, fc_threshold, P_VALUE, width, height)
        else:
            st.info(f"Значимых генов с p < {P_VALUE} and |FC| ≥ 1 нет.")


main()
