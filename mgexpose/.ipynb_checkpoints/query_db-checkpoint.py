import psycopg2
import pandas as pd
import json
import sys

genome2speci_pg3_file = '/home/grekova/workspace/promge_website/data/genome_id2speci.tsv'

df = pd.read_csv(genome2speci_pg3_file, sep="\t")

genome2speci_ids = dict(zip(df["sample_name"], df["cluster_name"]))

def get_gtdb_id(identifier):
    """Extract the GTDB ID from a genome identifier."""
    parts = identifier.split('_')
    return f"{parts[0]}_{parts[1]}"  # GCA_xxxxxx.x
    
def connect(params_dic):
    """ Connect to the PostgreSQL database server """
    conn = None
    try:
        # connect to the PostgreSQL server
        print('Connecting to the PostgreSQL database...')
        conn = psycopg2.connect(**params_dic)
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
        sys.exit(1) 
    print("Connection created successfully")
    return conn
    
def parse_mge_id(id):
    #'GCA_009102765.1_371601.SAMN11944272.WDCH01000111:267-2860'
    mge_dict = {'gtdb_id': '',
                'contig' : '',
                'start' : '',
                'end' : ''
               }
    id = id.replace(':', '_').split('_')
    mge_dict['gtdb_id'] = (id[0] + '_' + id[1]) #GCA_xxxxxx.x
    mge_dict['contig'] = (id[2])
    coordinates = [int(c) for c in id[3].split('-')]
    mge_dict['start'] = (coordinates[0])
    mge_dict['end'] = (coordinates[1])
    return mge_dict

def query_mge_annotations(conn):
    ''' Query the database to get all levels of taxonomy, mge_type and recombinases from the mge table.
    Args:
        In: conn: psycopg connection
        Out: result: df mge_id recombinase tax_domain tax_phylum tax_class tax_order tax_family tax_genus tax_species     

    '''
    cursor = conn.cursor()
    levels = ["clusters.{level}".format(level=level) for level in ['tax_domain', 'tax_phylum', 'tax_class', 
              'tax_order', 'tax_family', 'tax_genus', 
              'tax_species']]
    levels_str = ', '.join(levels)
    query = """
    SELECT contig || ':' || start_pos || '-' || end_pos AS contig_pos,
    {levels_str}
    FROM clusters AS clusters, pg3.mge AS mge
    WHERE clusters.id = mge.cluster_id;
    """.format(levels_str=levels_str)
    cursor.execute(query)

    result = cursor.fetchall()
    cursor.close()
    columns = ['contig_pos']
    columns.extend([l.replace("clusters.", "") for l in levels])
    result = pd.DataFrame(result, columns=columns)
    return result

def get_taxa(speci_lst, cursor, level=None):
    ''' Query the database to get taxonomy information
    Args:
        In: speci_lst (list): List of species names
            cursor: psycopg cursor object
            level (str, optional): Specific taxonomic level to query. If None, fetch all levels.
        Out: result: (DataFrame) containing taxonomy information
    '''
    levels = {
        "tax_domain", "tax_phylum", "tax_class", "tax_order", "tax_family", "tax_genus", "tax_species"
    }
    
    if level and level not in levels:
        raise ValueError(f"Invalid level: {level}. Choose from {levels} or None for full taxonomy.")
    
    levels_str = f"clusters.{level}" if level else ', '.join(f"clusters.{lvl}" for lvl in levels)
    
    specI_str = ', '.join(['%s'] * len(speci_lst))
    
    query = f"""
    SELECT cluster_name, {levels_str}
    FROM clusters AS clusters
    WHERE clusters.cluster_name IN ({specI_str});
    """
    
    cursor.execute(query, tuple(speci_lst))
    result = cursor.fetchall()
    
    columns = ['cluster_name'] + ([level] if level else list(levels))
    
    return pd.DataFrame(result, columns=columns) if result else pd.DataFrame(columns=columns)

def get_gtdb_taxa(sample_ids, db, cursor, level=None):
    levels = {
        "d", "p", "c", "o", "f", "g", "s"
    }
    if level and level not in levels:
        raise ValueError(f"Invalid level: {level}. Choose from {levels} or None for full taxonomy.")
    
    if db == "pg3":
        tax_table = "pg3.gtdb_r220"
        sample_table = "pg3.samples" 
        assembly = "genome_id"
        levels_str = f"t.{level}" if level else ', '.join(f"t.{lvl}" for lvl in levels)
    elif db == "spire":
        tax_table = "gtdb_r220"
        sample_table = "bins"
        assembly = "bin_id"
        levels_str = f"{tax_table}.{level}" if level else ', '.join(f"{tax_table}.{lvl}" for lvl in levels)
    else:
        raise ValueError(f"Invalid db specification: {db}. pg3 or spire are allowed.")
    
    
    sample_ids_str = ', '.join(['%s'] * len(sample_ids))
    
    if db == "pg3":   
        query = f"""
        SELECT sample_name, {levels_str}
        FROM {sample_table} AS s, {tax_table} AS t
        WHERE (s.sample_name IN ({sample_ids_str})) AND (s.id = t.sample_id);
        """
    elif db == "spire":
        query = f"""
        SELECT bin_name, {levels_str}
        FROM {sample_table} AS {sample_table}, {tax_table} AS {tax_table}
        WHERE ({sample_table}.bin_name IN ({sample_ids_str})) AND ({sample_table}.id = {tax_table}.bin_id);
        """
    else:
        raise ValueError(f"Invalid db specification: {db}. pg3 or spire are allowed.")
    
    cursor.execute(query, tuple(sample_ids))
    result = cursor.fetchall()
    
    columns = [assembly] + ([level] if level else list(levels))
    
    return pd.DataFrame(result, columns=columns) if result else pd.DataFrame(columns=columns)


def annotate_clustering_df(clustered_df, conn, level="tax_species"):
    ''' Query taxonomy information and update DataFrame '''
    
    print("Clustering df, nrows:", len(clustered_df))
    
    genome2speci = {id: genome2speci_ids[get_gtdb_id(id)] for id in clustered_df.member_seq_100 if 'GCA_' in id}
    clustered_df['speci'] = clustered_df['member_seq_100'].map(genome2speci)
    specIs = list(set(genome2speci.values()))
    print("# specI:", len(specIs))
    
    cursor = conn.cursor()
    
    if level == "full":
        result_df = get_taxa(specIs, cursor)
    else:
        result_df = get_taxa(specIs, cursor, level)
    
    print("Merging taxonomy")
    clustered_df = clustered_df.merge(result_df, how="inner", left_on="speci", right_on="cluster_name")
    
    cursor.close()
    return clustered_df


def get_speci_taxonomy_df(speci_lst, conn, level="tax_species"):
    ''' For each specI cluster get taxonomy information and return as taxa_df dataframe'''
    
    specIs = list(set(speci_lst)) # Ensure that it is a list
    print("# specI:", len(speci_lst))
    
    cursor = conn.cursor()
    
    if level == "full":
        taxa_df = get_taxa(speci_lst, cursor)
    else:
        taxa_df = get_taxa(speci_lst, cursor, level)
    
    cursor.close()
    return taxa_df


def get_gtdb_taxonomy_df(sample_ids, db, conn, level="tax_species"):
    ''' For each sample_id (bin_id or genome_id) get gtdb taxonomy information and return as taxa_df dataframe'''
    
    sample_ids= list(set(sample_ids)) # Ensure that it is a list
    print("# samples_ids:", len(sample_ids))

    cursor = conn.cursor()
    
    if level == "full":
        taxa_df = get_gtdb_taxa(sample_ids, db, cursor)
    else:
        taxa_df = get_gtdb_taxa(sample_ids, db, cursor, level)
    
    cursor.close()
    return taxa_df


def get_microontology(sample_names, conn):
    ''' 
    Query the database to get microontology information.
    
    Args:
        sample_names (list): List of sample names
        cursor (psycopg cursor): Active DB cursor
    
    Returns:
        pd.DataFrame: DataFrame with sample_name, sample_id, term_id, term, term_array
    '''

    if len(sample_names) == 0:
        return pd.DataFrame(columns=["sample_name", "study_id", "sample_id", "term"])
    cursor = conn.cursor()

    samples_str = ', '.join(['%s'] * len(sample_names))

    query = f"""
    SELECT 
        s.sample_name,
        s.study_id,
        mv.sample_id,
        mt.term
    FROM samples s
    JOIN microntology_v3 mv ON s.id = mv.sample_id
    JOIN LATERAL unnest(mv.microntology_terms) AS term_id ON TRUE
    JOIN microntology_terms mt ON mt.id = term_id
    WHERE s.sample_name IN ({samples_str});
    """

    cursor.execute(query, tuple(sample_names))
    result = cursor.fetchall()
    cursor.close()
    columns = ["sample_name", "study_id", "sample_id", "term"]
    
    return pd.DataFrame(result, columns=columns) if result else pd.DataFrame(columns=columns)




































