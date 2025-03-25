import streamlit as st
import pandas as pd



st.title("Gene Enrichment Analysis")

# Check if group datasets exist in session state
# if 'group_datasets' not in st.session_state or not st.session_state.group_datasets:
#     st.warning("Please load and prepare datasets on the Load Data page first!")
#     st.stop()

from enrichr_analyzer import EnrichrAnalyzer

# Ввод данных
gene_input = st.text_area("Enter gene symbols (one per line)", "BRCA1\nTP53\nEGFR\nMYC\nCDKN2A")
gene_list = [g.strip() for g in gene_input.split('\n') if g.strip()]

description = st.text_input("Analysis description")
library = st.selectbox("Select library", ["KEGG_2016", "GO_Biological_Process_2021"])

if len(gene_input) and st.button("Run Analysis"):
    analyzer = EnrichrAnalyzer()
    
    # Получаем результаты
    with st.spinner("Running enrichment analysis..."):
        results = analyzer.enrich(gene_list, description, library)
    
    # Показываем таблицу
    st.subheader("Enrichment Results")
    st.dataframe(results)
    
    st.markdown("""
    ### **Как создается сеть обогащения и что она показывает?**  

    Сеть обогащения (Enrichment Network) — это граф, который визуализирует связи между **генами** и **биологическими терминами** (например, путями KEGG или GO-терминами), полученными в результате анализа обогащения (например, через Enrichr).  

    ---

    ## **1. Структура сети**
    Граф состоит из **двух типов узлов** и **ребер** между ними:  

    ### **🔵 Узлы (Nodes)**
    - **🔶 Термины (Квадраты)**  
    - Это результаты обогащения (например, "Bladder cancer", "MicroRNAs in cancer").  
    - **Цвет** показывает значимость: чем **краснее**, тем **меньше p-value** (более значимый термин).  
    - **Размер** может отражать количество связанных генов или комбинированный балл (Combined Score).  

    - **🔵 Гены (Круги)**  
    - Это гены из вашего исходного списка, которые попали в данный термин.  
    - Обычно отображаются **синим** цветом.  

    ### **🔗 Ребра (Edges)**  
    - Показывают, **какие гены входят в какой термин**.  
    - Если ген связан с несколькими терминами, это видно по количеству соединений.  

    ---

    ## **2. Как строится граф?**  
    1. **Из результатов Enrichr берутся топовые термины** (например, 20 самых значимых).  
    2. **Для каждого термина извлекаются гены**, которые в него входят.  
    3. **Создается граф** (`networkx.Graph`), где:  
    - **Термины → квадраты**  
    - **Гены → круги**  
    - **Ребра → связи "термин-ген"**  
    4. **Раскраска и позиционирование:**  
    - Термины раскрашиваются по **p-value** (красный = значимый).  
    - Гены всегда **синие**.  
    - Позиции узлов рассчитываются алгоритмом (`spring_layout`), чтобы связанные элементы были ближе.  

    ---

    ## **3. Что можно увидеть в таком графе?**  
    ✅ **Какие гены входят в несколько путей** (например, TP53 может быть и в раке, и в апоптозе).  
    ✅ **Какие термины самые значимые** (чем краснее, тем важнее).  
    ✅ **Есть ли кластеры** (группы терминов, связанных через общие гены).  

    ---

    ## **4. Пример интерпретации**  
    Допустим, у вас есть граф, где:  
    - **Термин "Bladder cancer"** (красный квадрат) связан с генами **TP53, EGFR**.  
    - **Термин "MicroRNAs in cancer"** (оранжевый квадрат) связан с **TP53, MYC**.  

    **Вывод:**  
    - **TP53** участвует в обоих процессах → возможно, ключевой ген.  
    - **"Bladder cancer"** более значим (краснее), чем **"MicroRNAs in cancer"**.  

    ---

    ### **Итог**  
    Сеть обогащения помогает **визуально** оценить:  
    1. **Какие гены наиболее важны** (много связей).  
    2. **Какие биологические процессы/пути наиболее релевантны** (цвет и размер узлов).  
    3. **Есть ли общие паттерны** (кластеры терминов).  

    Это гораздо удобнее, чем просто таблица с p-value! 🚀
    """)
    
    # Показываем сеть
    st.subheader("Enrichment Network")
    analyzer.plot_enrichment_network(results)
