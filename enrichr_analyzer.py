import json
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import streamlit as st
import numpy as np


class EnrichrAnalyzer:
    """
    Расширенный класс для анализа обогащения генов с визуализацией сети
    """

    BASE_URL = 'http://amp.pharm.mssm.edu/Enrichr/'

    def __init__(self):
        pass

    def _add_gene_list(self, gene_list: list, description: str) -> dict:
        """Добавляет список генов в Enrichr"""
        genes_str = '\n'.join(gene_list)
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }
        response = requests.post(self.BASE_URL + 'addList', files=payload)
        if not response.ok:
            raise Exception(f'Error adding gene list: {response.status_code}')
        return json.loads(response.text)

    def _get_enrichment_results(self, user_list_id: str, library: str) -> dict:
        """Получает результаты обогащения"""
        query = f'enrich?userListId={user_list_id}&backgroundType={library}'
        response = requests.get(self.BASE_URL + query)
        if not response.ok:
            raise Exception(
                f'Error getting enrichment results: {response.status_code}')
        return json.loads(response.text)

    def enrich(self, gene_list: list, description: str, library: str,
               top_terms: int = 20) -> pd.DataFrame:
        """
        Выполняет анализ обогащения с ограничением по количеству терминов

        Args:
            top_terms: количество топовых терминов для возврата
        """
        job_id = self._add_gene_list(gene_list, description)
        results = self._get_enrichment_results(job_id['userListId'], library)
        df = self._parse_enrichment_results(results, library)
        return df.head(top_terms)

    def _parse_enrichment_results(self, results: dict, library: str) -> pd.DataFrame:
        """Парсит JSON-результаты в DataFrame"""
        if library not in results:
            raise ValueError(f"Library {library} not found in results")

        records = []
        for item in results[library]:
            records.append({
                'Rank': item[0],
                'Term': item[1],
                'P-value': item[2],
                'Z-score': item[3],
                'Combined Score': item[4],
                'Genes': item[5],
                'Adjusted P-value': item[6],
                'Old P-value': item[7],
                'Old Adjusted P-value': item[8]
            })

        df = pd.DataFrame(records)
        df['-log10(P-value)'] = -np.log10(df['P-value'])
        return df.sort_values('P-value')

    def plot_enrichment_network(self, df: pd.DataFrame, figsize=(10, 10)):
        """
        Строит сеть обогащения и отображает в Streamlit

        Args:
            df: DataFrame с результатами обогащения
            figsize: размер фигуры
        """
        G = nx.Graph()

        # Добавляем узлы для терминов и генов
        for _, row in df.iterrows():
            term = row['Term']
            G.add_node(term, type='term', pvalue=row['P-value'])

            for gene in row['Genes']:
                G.add_node(gene, type='gene')
                G.add_edge(term, gene)

        # Раскраска узлов
        term_nodes = [n for n, attrs in G.nodes(
            data=True) if attrs['type'] == 'term']
        gene_nodes = [n for n, attrs in G.nodes(
            data=True) if attrs['type'] == 'gene']

        # Позиционирование узлов
        pos = nx.spring_layout(G, k=0.3, iterations=50)

        # Создаем фигуру
        fig, ax = plt.subplots(figsize=figsize)

        # Рисуем термины (круги)
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=term_nodes,
            node_color=[G.nodes[n]['pvalue'] for n in term_nodes],
            node_size=800,
            cmap=plt.cm.Reds,
            alpha=0.8,
            node_shape="s"
        )

        # Рисуем гены (квадраты)
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=gene_nodes,
            node_color='skyblue',
            node_size=400,
            alpha=0.8,
            node_shape="o"
        )

        # Рисуем ребра
        nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.2)

        # Подписи узлов
        nx.draw_networkx_labels(
            G, pos,
            labels={n: n for n in G.nodes()},
            font_size=8,
            font_family='sans-serif'
        )

        # Легенда и цветовая шкала
        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.Reds,
            norm=plt.Normalize(vmin=df['P-value'].min(), vmax=df['P-value'].max()))
        sm.set_array([])
        plt.colorbar(
            sm,
            ax=ax,  # Указываем конкретный Axes
            label='P-value',
            shrink=0.5,
            location='right'
        )

        plt.title("Enrichment Network", fontsize=12)
        plt.axis('off')

        # Отображаем в Streamlit
        st.pyplot(plt.gcf())
        plt.clf()
