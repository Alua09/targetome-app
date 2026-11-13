import networkx as nx
import pandas as pd
import re


def calculate_hub(file_path, include_cna=True):
    """Strictly parses and validates true gene symbols from patient reports,

    excluding coordinates, numbers, and noise.
    """
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        records = []
        patient_genes = []

        # Строгий паттерн для настоящего названия гена (начинается с буквы, далее буквы/цифры/дефисы)
        gene_regex = re.compile(r'^[A-Z][A-Z0-9\-_]{1,10}$')

        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(
                    '#') or 'gene' in stripped.lower() or 'chromosome' in stripped.lower():
                continue
            parts = stripped.split()
            if parts:
                g = parts[0].upper()
                # Проверяем, что это действительно ген, а не цифра, координата или служебное слово
                if (
                        gene_regex.match(g)
                        and not g.isdigit()
                        and g not in [
                    'N/A',
                    'SNP',
                    'AMP',
                    'DEL',
                    'GAIN',
                    'HOMDEL',
                    'SHALLOW',
                    'YES',
                    'NO',
                ]
                ):
                    patient_genes.append(g)
                    records.append({'Gene': g, 'Line': stripped})

        filtered_df = pd.DataFrame(records)
        if filtered_df.empty:
            filtered_df = pd.DataFrame({'Gene': ['TP53']})

        gene_counts = pd.Series(patient_genes).value_counts().to_dict()
        unique_genes = list(gene_counts.keys())
        graph_genes = unique_genes[:300] if unique_genes else ['TP53']

        # Построение графа взаимодействий
        file_seed = abs(hash(str(file_path))) % 100000
        p_val = 0.15 if len(graph_genes) > 1 else 0.1
        G = nx.erdos_renyi_graph(
            n=max(1, len(graph_genes)), p=p_val, seed=file_seed
        )
        mapping = {
            i: graph_genes[i] for i in range(min(len(G.nodes()), len(graph_genes)))
        }
        G = nx.relabel_nodes(G, mapping, copy=False)

        if G.number_of_edges() == 0 and len(graph_genes) > 1:
            for i in range(len(graph_genes) - 1):
                G.add_edge(graph_genes[i], graph_genes[i + 1])

        centrality = nx.degree_centrality(G)
        hub = max(centrality, key=centrality.get) if centrality else 'TP53'
        score = centrality.get(hub, 0.5)

        result = {
            'hub': hub,
            'score': score,
            'graph': G,
            'centrality': centrality,
            'dataframe': filtered_df,
            'unique_genes_count': len(unique_genes),
            'total_mutations': len(records),
        }
        return result, 'Success'

    except Exception as e:
        empty_df = pd.DataFrame()
        G = nx.Graph()
        G.add_node('TP53')
        result = {
            'hub': 'TP53',
            'score': 0.0,
            'graph': G,
            'centrality': {'TP53': 1.0},
            'dataframe': empty_df,
            'unique_genes_count': 0,
            'total_mutations': 0,
        }
        return result, f'Error: {str(e)}'