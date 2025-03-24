import streamlit as st

# Настройка страницы
st.set_page_config(
    page_title="Differential Gene Expression Analysis",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# Главный заголовок
st.title("🧬 Differential Gene Expression Analysis")
st.markdown("---")

# Описание программы
st.markdown("""
## 🌟 О программе

Этот инструмент предназначен для анализа дифференциальной экспрессии генов (DGE) 
из датасетов [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).

Программа позволяет:
- Загружать и предобрабатывать данные RNA-seq или микрочипов
- Проводить анализ качества данных
- Выявлять дифференциально экспрессируемые гены
- Визуализировать результаты анализа
""")

st.markdown("---")

# Возможности
st.markdown("""
## 🔍 Основные возможности

### 📊 Загрузка данных
- Поддержка форматов GEO (Series Matrix, SOFT)
- Автоматическое извлечение метаданных
- Нормализация сырых данных

### 🔬 Анализ дифференциальной экспрессии
- Методы: DESeq2, edgeR, limma-voom (для RNA-seq)
- Коррекция множественного тестирования (FDR)
- Анализ обогащения (GO, KEGG)

### 📈 Визуализация
- Heatmaps экспрессии генов
- Volcano plots
- PCA plots для оценки группировки образцов
""")

st.markdown("---")

# Как использовать
st.markdown("""
## 🚀 Как использовать

1. **Загрузите данные** через интерфейс или укажите GEO accession number
2. **Настройте параметры** анализа
3. **Запустите анализ** и дождитесь результатов
4. **Исследуйте результаты** через интерактивные графики
5. **Экспортируйте** таблицы и изображения

Для начала работы перейдите в соответствующую вкладку в боковом меню.
""")

st.markdown("---")

# Контактная информация
st.markdown("""
## 📧 Контакты

По вопросам использования программы обращайтесь:

- Email: bioinfo@example.com
- GitHub: [github.com/your_repo](https://github.com/your_repo)
""")

# Боковая панель
with st.sidebar:
    st.title("🔬 DGE Analysis")
    st.markdown("""
    ### Навигация
    - **Главная** (текущая страница)
    - **Загрузка данных**
    - **Анализ DGE**
    - **Визуализация**
    - **Экспорт результатов**
    """)
    
    st.markdown("---")
    st.markdown("### Информация")
    st.info("""
    Версия: 1.0.0  
    Последнее обновление: 2023-11-15  
    Лицензия: MIT
    """)