import streamlit as st
import os
import pandas as pd
import subprocess
from pathlib import Path
from datetime import datetime

# ========== Config ==========
DATA_DIR = "data"
Path(DATA_DIR).mkdir(exist_ok=True)
st.title("Загрузка GEO-файлов")

# ========== Utility Functions ==========


@st.cache_data(ttl=5)
def get_files(extension='.gz'):
    files = []
    for item in os.listdir(DATA_DIR):
        path = os.path.join(DATA_DIR, item)
        if os.path.isfile(path) and item.endswith(extension):
            files.append({
                "Имя файла": item,
                "Размер (KB)": round(os.path.getsize(path) / 1024, 2),
                "Изменено": datetime.fromtimestamp(os.path.getmtime(path))
            })
    if not len(files):
        return pd.DataFrame()
    return pd.DataFrame(files).sort_values("Изменено", ascending=False)


def extract_name_from_path(path):
    return path.split("\\")[1].split("_")[0]


def run_r_script(input_path, name):
    return subprocess.run(
        ["Rscript", "process_geo.R", input_path, name],
        capture_output=True,
        text=True,
        check=True
    )


def read_expression_data(name):
    path = f"{DATA_DIR}//{name}//{name}_expr.csv"
    return pd.read_csv(path, index_col=0).rename_axis('id')


def read_phenotype_data(name):
    path = f"{DATA_DIR}//{name}//{name}_phen.csv"
    return pd.read_csv(path, index_col=0).rename_axis('id')


def read_gene_list(txt_path):
    with open(txt_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# ========== Display Functions ==========


def display_dataframe(title, df):
    st.subheader(title)
    st.dataframe(df)


def file_selector(title, df, key):
    st.subheader(title)
    return st.dataframe(
        df,
        key=key,
        on_select="rerun",
        selection_mode="single-row",
        use_container_width=True
    )


def display_debug_info(gz_selection):
    with st.expander("Debug Info"):
        st.write("GZ Selection object:", gz_selection)
        st.write("TXT Files DataFrame:", get_files('.txt'))

# ========== Main Logic ==========


def main():
    files_df = get_files('.gz')
    if files_df.empty:
        st.warning("GZ файлы не найдены в папке 'data'")
        return

    gz_selection = file_selector(
        "Доступные GEO датасеты", files_df, "gz_file_selector")

    if not gz_selection['selection']['rows']:
        st.info("👆 Выберите GEO-датасет для обработки")
        return

    selected_row = gz_selection['selection']['rows'][0]
    selected_file = files_df.iloc[selected_row]["Имя файла"]
    input_path = os.path.join(DATA_DIR, selected_file)
    name = extract_name_from_path(input_path)

    with st.spinner(f"Processing {selected_file}..."):
        try:
            output_expr_path = f"{DATA_DIR}//{name}//{name}_expr.csv"
            if not os.path.exists(output_expr_path):
                run_r_script(input_path, name)
                st.success(f"Successfully processed {selected_file}!")

            if os.path.exists(output_expr_path):
                expr_df = read_expression_data(name)
                phen_df = read_phenotype_data(name)

                display_dataframe("Матрица экспрессии", expr_df)
                display_dataframe("Данные фенотипов", phen_df)

                handle_gene_list_and_filtering(expr_df, phen_df)

        except subprocess.CalledProcessError as e:
            st.error(f"Processing failed: {e.stderr}")

    display_debug_info(gz_selection)

# ========== Gene List Selection and Filtering ==========


def handle_gene_list_and_filtering(expr_df, phen_df):
    st.markdown("---")
    txt_files_df = get_files('.txt')

    if txt_files_df.empty:
        st.warning("Не было найдено TXT файлов в папке 'data'")
        return

    txt_selection = file_selector(
        "Выбор списка генов для фильтрации", txt_files_df, "txt_file_selector")
    select_all = st.checkbox(
        "Выбрать все гены", value=False, key="select_all_genes")

    filtered_expr = expr_df
    if not select_all and txt_selection['selection']['rows']:
        row = txt_selection['selection']['rows'][0]
        txt_file = txt_files_df.iloc[row]["Имя файла"]
        txt_path = os.path.join(DATA_DIR, txt_file)
        gene_list = read_gene_list(txt_path)

        filtered_expr = expr_df[expr_df.index.isin(gene_list)]
        st.subheader(f"Отфильтрованная матрица экспресси (по {txt_file})")
        st.write(f"Генов в списке: {len(gene_list)}")
        st.write(f"Генов, найденных в матрице: {len(filtered_expr)}")
        st.dataframe(filtered_expr)

    if not filtered_expr.empty:
        handle_phenotype_filtering(filtered_expr, phen_df)

# ========== Phenotype Filtering and Group Creation ==========


def handle_phenotype_filtering(filtered_expr, phen_df):
    st.markdown("---")
    st.subheader("Группировка по фенотипу")

    phen_columns = [
        col for col in phen_df.columns if col != phen_df.index.name]
    if not phen_columns:
        st.warning("Не найдено колонок с фенотипами для группировки")
        return

    selected_col = st.selectbox(
        "Выберите колонку фенотипа для группировки:", phen_columns, key="phen_column_selector")
    values = phen_df[selected_col].unique()
    selected_values = st.multiselect(
        f"Выберите значения из {selected_col}, которые нужно оставить:", values, default=values, key="phen_value_selector")

    if not selected_values:
        return

    filtered_samples = phen_df[phen_df[selected_col].isin(
        selected_values)].index.tolist()
    final_filtered = filtered_expr[filtered_samples]

    display_dataframe("Финальная матрица экспрессии", final_filtered)

    if st.button("Создать датасеты по группам", key="create_group_datasets"):
        create_group_datasets(filtered_expr, phen_df, selected_col)

    if 'group_datasets' in st.session_state and st.session_state.group_datasets:
        display_group_datasets()

# ========== Group Dataset Logic ==========


def create_group_datasets(expr_df, phen_df, phen_column):
    groups = phen_df[phen_column].unique()
    group_datasets = {}
    avg_group_datasets = {}

    for group in groups:
        sample_ids = phen_df[phen_df[phen_column] == group].index.tolist()
        group_data = expr_df[sample_ids]
        group_datasets[group] = group_data
        avg_group_datasets[group] = pd.DataFrame(
            {'Expression': group_data.mean(axis=1)})

    combined_df = pd.concat([
        df.reset_index()
        .assign(Group=name)
        .rename(columns={'index': 'id'})
        for name, df in avg_group_datasets.items()
    ])
    combined_df = combined_df.pivot_table(
        index='Group', columns='id', values='Expression')

    st.session_state.group_datasets = {
        'groups': list(group_datasets.keys()),
        'datasets': group_datasets,
        'avg': combined_df
    }
    st.success(f"Создано {len(groups)} датасетов по группам!")


def display_group_datasets():
    st.markdown("---")
    st.subheader("Датасеты по группам")
    groups = st.session_state.group_datasets['groups']
    datasets = st.session_state.group_datasets['datasets']
    avg_df = st.session_state.group_datasets['avg']

    tabs = st.tabs([f"Группа {group}" for group in groups])
    for tab, group in zip(tabs, groups):
        with tab:
            st.write(f"Образцов в группе: {len(datasets[group].columns)}")
            st.dataframe(datasets[group])

    st.markdown("---")
    st.subheader("Датасет со средними значениями экспрессии по группам")
    st.dataframe(avg_df)


# Run the app
if __name__ == "__main__":
    main()
