import os
import pandas as pd


def collect_functions(WORKDIR, anno_path):
    '''
    Collects functions from BAKTA or PROKKA annotation files
    WORKDIR: str, path to the working directory
    anno_path: str, path to the annotation files (e.g. processed/bakta)
    Retuns a pandas DataFrame with columns: locus, genome_id, product, go annotation
    '''
    
    from glob import glob
    anno_folders = glob(os.path.join(WORKDIR, anno_path, "*"))
    all_data = pd.DataFrame()
    for anno_folder in anno_folders:
        genome_id = anno_folder.split(os.sep)[-1]
        functions = os.path.join(anno_folder, f'{genome_id}.tsv')
        data = pd.read_csv(functions, sep='\t', skiprows=5)
        data['genome_id'] = genome_id
        data = data[['Product', 'DbXrefs', 'Locus Tag', 'genome_id']]
        all_data = pd.concat([all_data, data])
    all_data = all_data.rename(columns={'Product':'product','DbXrefs':'go', 'Locus Tag': 'locus'})
    print(f'loaded {len(all_data)} functions from {len(anno_folders)} BAKTA tsv files')
    return all_data[['locus','genome_id','product','go']].drop_duplicates()


def get_pg_to_locus_map(WORKDIR, SPECIES):
    '''
    This function reads the allele names file and returns a dataframe with the mapping of pg cluster to locus in each strain
    '''
    fin = open(os.path.join(WORKDIR, f'processed/cd-hit-results/{SPECIES}_allele_names.tsv'))
    out = []
    for line in fin:
        allele = line.strip().split('\t')[0]
        cluster = allele.split('A')[0]
        loci = line.strip().split('\t')[1:]
        for locus in loci:
            out.append({'cluster': cluster, 'gene_id': locus})
    out = pd.DataFrame(out)
    out.drop_duplicates(inplace=True)
    return pd.DataFrame(out)


def calc_enrichment(L, cluster_to_go_functions, function, functions2genes, phylon, phylon_contribution_cutoff=0):
    '''
    calculates the enrichment of a function in a phylon
    input: 
        L = L matrix from NMF
        cluster_to_go_functions dataframe from get_cluster_to_go_functions function
        function = function of interest
        phylon = phylon of interest
        cutoff = cutoff for the phylon contribution
        
    output: p-value of the enrichment

    '''
    from scipy import stats

    clusters_in_phylon = L[L[phylon]>phylon_contribution_cutoff][phylon].index # TODO consider other ways to generate this cutoff (e.g. quantile, mass %, etc)
    
    clusters_w_function= cluster_to_go_functions[cluster_to_go_functions['go'] == function]['cluster'].values.tolist()
    len(clusters_w_function)
    overlap = set(clusters_in_phylon).intersection(set(clusters_w_function))
    common_proteins = '; '.join(functions2genes[functions2genes['cluster'].isin(overlap)]['product'].values.tolist())
    N = len(cluster_to_go_functions['cluster'].drop_duplicates()) # Total number of clusters
    K = len(clusters_w_function) # clusters with function in poi
    n = len(clusters_in_phylon) # Number of clusters in the phylon
    k = len(overlap) # Number of genes in the phylon with function in poi
    p_value = stats.hypergeom.sf(k - 1, N, K, n)
    
    return {'function': function, 'N':N, 'K':K, 'n':n, 'k':k, 'p_value':p_value, 'products':common_proteins}


def calc_all_phylon_go_enrichments(L, all_functions, cluster_to_go_functions, go_functions, phylon_contribution_cutoff=0):
    '''
    calculates the enrichment of all functions in all phylons
    input: 
        L = L matrix from NMF
        cluster_to_go_functions dataframe from get_cluster_to_go_functions function
        go_functions = go functions dataframe
        cutoff = cutoff for the phylon contribution
    '''
    import numpy as np
    out = []
    for phylon in L.columns:
        for function in go_functions.index:
            res = calc_enrichment(L, cluster_to_go_functions, function, all_functions, phylon, phylon_contribution_cutoff=phylon_contribution_cutoff)
            # print(function, res['p_value'])
            out.append({
                'phylon': phylon,
                'function': function,
                'p_value': res['p_value'],
                'genes in phylon': res['n'],
                'genes w function': res['K'],
                'overlap': res['k'],
                'products': res['products']
            })
    out = pd.DataFrame(out).sort_values('p_value')
    out['logp'] = -1 * out['p_value'].apply(np.log10)

    #TODO: add FDR

    return out

def get_go_mapping():
    '''
    returns a dataframe of go terms
    loads from the data/go_terms.csv file - need to keep this updated
    '''
    
    # go_mapping = os.path.join('../pyphylon/data', 'go_terms.csv')
    go_mapping = os.path.join(os.sep.join(__file__.split(os.sep)[:-2]), 'pyphylon', 'data', 'go_terms.csv')
    print(f'loaded go terms from {go_mapping}')
    go_mapping = pd.read_csv(go_mapping, index_col=0)
    return go_mapping


def explode_go_annos(functions):
    '''
    generates a dataframe with the cluster, phylon, and function
    input: L matrix from NMF and functions dataframe and a function dataframe
    
    output: dataframe with cluster, phylon, and function
    '''
    acc_functions = functions # don't filter by accessory genome...

    cluster_to_go_functions = acc_functions[['cluster','go']].copy()
    cluster_to_go_functions['go'] = cluster_to_go_functions['go'].astype(object).str.split(', ')
    cluster_to_go_functions = cluster_to_go_functions.explode('go')
    cluster_to_go_functions
    
    return cluster_to_go_functions.reset_index(drop=True)



def gen_phylon_wordcloud(L, functions, phylon, cutoff=0,  save=False, filename='phylon_wordcloud.png'):
    '''
    generates a wordcloud of the functions in a phylon
    input: 
        L = L matrix from NMF
        functions = functions dataframe
        phylon = phylon of interest
        cutoff = cutoff for the phylon contribution
    '''
    from wordcloud import WordCloud
    import matplotlib.pyplot as plt

    coi = L[L[phylon]>cutoff][phylon].index
    
    cluster_functions = functions[functions['cluster'].isin(coi)]['product']
    # print(cluster_functions.tolist())

    # Combine all terms into a single string
    terms_string = ' '.join(cluster_functions).lower() \
    .replace('protein','') \
    .replace('hypothetical', '') \
    .replace('uncharacterized','') \
    .replace('putative','') \
    .replace('domain','') \
    .replace('containing','') \
    .replace('family','') \
    .replace('like','') \
    .replace('protein','') \
    .replace('related','') \
    .replace('conserved','') \
    .replace('possible','') \
    .replace('unknown','') \
    .replace('transcriptional regulator','') \
    .replace('subunit','') \
    .replace('type','') \
    .replace('system','')
        
    

    # Generate the word cloud
    wordcloud = WordCloud(width=1200, height=1200, background_color='white').generate(terms_string)

    # Display the word cloud
    plt.figure(figsize=(15, 15))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.title(phylon)
    if save:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()
        plt.close()


def load_eggnog_annotations(data_dir, species):
    """Load eggNOG-mapper annotations and map allele IDs to cluster IDs.

    Parameters
    ----------
    data_dir : str
        Path to data directory (e.g. "data")
    species : str
        Species prefix (e.g. "CJejuni")

    Returns
    -------
    pd.DataFrame
        Indexed by cluster ID with columns: cog_category, kegg_ko,
        kegg_pathway, kegg_module, eggnog_desc
    """
    import re

    path = os.path.join(data_dir, "processed", "eggnog", f"{species}.emapper.annotations")
    df = pd.read_csv(path, sep="\t", comment="#", header=None)

    # eggNOG-mapper v2 output columns (0-indexed):
    # 0=query, 4=evalue, 5=score, 6=eggNOG_OGs, 7=max_annot_lvl,
    # 8=COG_category, 9=Description, 10=Preferred_name,
    # 11=GOs, 12=EC, 13=KEGG_ko, 14=KEGG_Pathway, 15=KEGG_Module,
    # 16=KEGG_Reaction, 17=KEGG_rclass, 18=BRITE, 19=KEGG_TC,
    # 20=CAZy, 21=BiGG_Reaction, 22=PFAMs
    col_map = {0: "query", 8: "cog_category", 9: "eggnog_desc",
               13: "kegg_ko", 14: "kegg_pathway", 15: "kegg_module"}
    df = df.rename(columns=col_map)[[*col_map.values()]]

    # Strip allele suffix: CJejuni_C1093A0 -> CJejuni_C1093
    df["cluster"] = df["query"].apply(lambda x: re.sub(r"A\d+$", "", x))
    df = df.drop(columns=["query"]).set_index("cluster")

    # Replace eggNOG-mapper's "-" placeholder with empty string
    df = df.replace("-", "")

    return df


def load_cogclassifier_results(data_dir, species):
    """Load COGclassifier results and map allele IDs to cluster IDs.

    Parameters
    ----------
    data_dir : str
        Path to data directory (e.g. "data")
    species : str
        Species prefix (e.g. "CJejuni")

    Returns
    -------
    pd.DataFrame
        Indexed by cluster ID with columns: cog_id, cog_category,
        cog_description
    """
    import re

    result_dir = os.path.join(data_dir, "processed", "cogclassifier")
    # COGclassifier outputs classifier_result.tsv in the output directory
    path = os.path.join(result_dir, "classifier_result.tsv")
    df = pd.read_csv(path, sep="\t")

    # Columns: QUERY_ID, COG_ID, CDD_ID, EVALUE, IDENTITY, GENE_NAME,
    #          COG_NAME, COG_LETTER, COG_DESCRIPTION
    df = df.rename(columns={
        "QUERY_ID": "query",
        "COG_ID": "cog_id",
        "COG_LETTER": "cog_category",
        "COG_DESCRIPTION": "cog_description",
    })

    # Strip allele suffix
    df["cluster"] = df["query"].apply(lambda x: re.sub(r"A\d+$", "", x))
    df = df[["cluster", "cog_id", "cog_category", "cog_description"]].set_index("cluster")

    return df
