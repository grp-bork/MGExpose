#!/usr/bin/env python
# pylint: disable=R0912,R0914,R0915
''' Collection of functions to collect MGEs from batches of files'''
import os
from collections import Counter, defaultdict
import json

from dask.distributed import Client, progress, WorkerPlugin
import dask.bag as db
from dask.bag import from_delayed
from dask.delayed import delayed
import dask

from .gffio import read_mge_genomic_islands_gff

from .base_logger import logger
import traceback


# Create a dictionary: genome: list of MGE IDs. This is needed to filter out only relevant MGEs per genome.The input is a corresponding list (same order) of genomes and MGE IDs.
def get_genome2mges(genomes, mges):
    genome2mge_id = {}
    for id, genome_id in zip(mges, genomes):
        if genome_id not in genome2mge_id:
            genome2mge_id[genome_id] = []  # Initialize the list
        genome2mge_id[genome_id].append(id)
    return genome2mge_id

# Helper function to extract genome ID/bin ID from file path
def get_genome_id_from_path(path):
    """
    Extract genome ID from a file path.
    """
    genome_id = None
    try:
        genome_id = path.split("/")[-2]
    except Exception as e:
        logger.error(f"Error extracting genome/bin ID{path}: {e}")
    return genome_id

def collect_batch_mges(gff_paths, i, relevant_ids=None):
    """
    Collect MGEs from a batch of GFF files.
    
    Parameters:
    - gff_paths: List of GFF file paths.
    - i: Index of the batch.
    - relevant_ids: Optional dictionary of relevant MGE IDs per genome ID.

    Returns:
    - List of MGE islands for all files in the batch.
    """
    islands = []
    for gff_path in gff_paths:
        genome_id = get_genome_id_from_path(gff_path)
        #logger.info(f"Processing genome: {genome_id}")

        try:
            if relevant_ids:
                relevant_mges = list(read_mge_genomic_islands_gff(gff_path, relevant_ids[genome_id]))
            else:
                relevant_mges = list(read_mge_genomic_islands_gff(gff_path))
            
            islands.extend(relevant_mges)

        except Exception as e:
            logger.error(f"Error processing {gff_path}: {e}")
            logger.error(traceback.format_exc())  # Full traceback

    logger.info(f"Batch {i} completed, MGE islands found: {len(islands)}")
    return islands


def apply_per_batch(islands, funcs):
    """
    Calculate statistics for a batch of MGEs using a list of functions.

    Parameters:
    - islands: List of MGE objects in the batch.
    - funcs: List of functions to be applied to the batch.

    Returns:
    - Dictionary with function names as keys and their results as values.
    """
    results = {}
    for func in funcs:
        func_name = func.__name__  # Get the function's name
        results[func_name] = func(islands)  # Apply the function and store the result
    return results

def apply_one_per_batch(islands, func):
    """
    Apply a single function to a batch of MGE objects.

    Parameters:
    - islands: List of MGE objects in the batch.
    - func: A single function to be applied to the batch.

    Returns:
    - The result of applying the function to the islands.
    """
    try:
        return func(islands)
    except Exception as e:
        raise RuntimeError(f"Error applying function '{func.__name__}' to batch: {e}")


def write_batch_json(batch_count, i, dir, base_filename):
    """
    Saves batch statistics to a JSON file using a Dask delayed function. This allows the function to be part of a larger
    Dask computation graph, potentially executed in parallel.

    Parameters:
    - batch_count (dict): Batch statistics, typically a dictionary with `Counter` values.
    - i (int): Index of the current batch. This is used to generate a unique filename for each batch.
    - dir (str): Path to the directory where the batch JSON files will be stored.
    - base_filename (str): Base name for the JSON files. The batch index will be appended to this base name to create the full filename.

    Ensures that the specified directory exists before writing the file. If it does not exist, the directory will be created.
    
    Returns:
    - A Dask delayed object which, when executed, will write the batch statistics to the specified file path in JSON format.
    """
    # Ensure directory exists
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Delayed function to write JSON
    def delayed_write(path, data):
        with open(path, 'w') as file:
            json.dump(data, file, indent=4)
    
    # Construct the full path with batch number
    path = os.path.join(dir, f"{base_filename}_{i}.json")
    return delayed(delayed_write)(path, batch_count)


def write_batch_tsv(tsv_string, i, dir, base_filename):
    """
    Saves a TSV-formatted string to a file as part of a Dask computation graph.

    Parameters:
    - tsv_string (str): TSV-formatted data returned from a processing function.
    - i (int): Index of the current batch, used to make unique filenames.
    - dir (str): Directory where TSV files will be saved.
    - base_filename (str): Base name for the output files.

    Returns:
    - A Dask delayed object that writes the TSV to disk when executed.
    """
    # Ensure the output directory exists
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Define the write function
    def delayed_write(path, content):
        with open(path, 'w') as file:
            file.write(content)

    # Create the output file path
    path = os.path.join(dir, f"{base_filename}_{i}.tsv")

    # Return the delayed write operation
    return delayed(delayed_write)(path, tsv_string)
    

def aggregate_attr(batches):
    """
    Aggregate string attributes across all batches.

    Parameters:
    - batches: List of batch statistics (dictionaries with str values e.g. cluster: COG category).

    {'90371.SAMN11043730.GCA_007435405_02914': 'S', '90371.SAMN14863315.GCA_013264555_00909': 'S', '28150.SAMN09228819.GCA_007140965_00837': 'K', '28901.SAMN13391507.GCA_011477875_00875': 'no_cog_fcat', '1967657.SAMN09203654.GCA_010924635_01295': 'S', '28901.SAMN13057743.GCA_009231785_02754': 'no_cog_fcat', '90371.SAMN11355433.GCA_007687065_04488': 'S', '28901.SAMN06645026.GCA_009179045_04086': 'no_cog_fcat', '28901.SAMN15658059.GCA_013797405_02028': 'S', '796732.SAMN01805325.GCA_000272735_04407': 'no_cog_fcat', '28901.SAMN12571445.GCA_010939435_00055': 'S', '28901.SAMN13747386.GCA_010741235_03029': 'no_cog_fcat', '115981.SAMN14080650.GCA_011486465_01735': 'S', '28901.SAMN14341880.GCA_011465135_00878': 'H', '1151002.SAMN09403228.GCA_004177825_01940': 'no_cog_fcat', '1029990.SAMN02415182.GCA_000484355_00589': 'S', '28901.SAMN12287151.GCA_007468615_01067': 'S', '28901.SAMN13057273.GCA_009230165_00580': 'no_cog_fcat', '611.SAMN21335643.GCA_019899165_00890': 'EH', '28901.SAMN10095790.GCA_005443695_02309': 'no_cog_fcat', '340190.SAMN15147492.GCA_013661605_01048': 'S', '224729.SAMN19336595.GCA_018502705_02901': 'no_cog_fcat', '28901.SAMN16355443.GCA_015155595_01828': 'no_cog_fcat', '59201.SAMN10093771.GCA_007777665_02600': 'no_cog_fcat', '59201.SAMN17835677.GCA_017072195_01803': 'S'}
    
2025-01-31 18:35:27,669 - {'611.SAMN17086052.GCA_016740915_04210': 'K', '611.SAMN07152477.GCA_007233055_03707': 'S', '1173835.SAMN01088029.GCA_000962395_02473': 'L', '1620419.SAMN03894126.GCA_001241425_04679': 'T', '568709.SAMEA2272227.GCA_000493535_02720': 'no_cog_fcat', '28901.SAMN18448990.GCA_017574325_01595': 'L', '90371.SAMEA6057931.GCA_016228905_04588': 'L', '28901.SAMN10177571.GCA_005772365_02802': 'N', '28901.SAMN14050865.GCA_011246635_04142': 'G', '90371.SAMN09387768.GCA_007158225_04690': 'N', '28144.SAMN07734943.GCA_003548115_02128': 'no_cog_fcat', '90105.SAMN09474912.GCA_004184575_03995': 'G', '59201.SAMN10756627.GCA_007583145_03925': 'S', '90371.SAMN03169328.GCA_008018515_03842': 'K', '1620419.SAMN04255380.GCA_010457935_04445': 'no_cog_fcat', '28901.SAMN16124589.GCA_014542005_01530': 'G', '28901.SAMN17005521.GCA_015838815_01302': 'no_cog_fcat', '28901.SAMN19285790.GCA_018468945_04349': 'G', '28901.SAMN10425133.GCA_010255835_04194': 'S', '28901.SAMN12344366.GCA_007726245_00622': 'S', '28901.SAMN12107692.GCA_006482085_01260': 'S', '440524.SAMN02867573.GCA_010663445_04243': 'M', '28901.SAMN20181473.GCA_020012815_04294': 'L', '28901.SAMEA6514879.GCA_011786425_00695': 'G', '90371.SAMN07279560.GCA_002260995_02309': 'K', '90371.SAMN19798225.GCA_018997055_00491': 'L', '28901.SAMN12823265.GCA_008717395_04569': 'V', '1173837.SAMN01088030.GCA_000962405_02064': 'L', '399584.SAMN13050934.GCA_009225065_04346': 'G', '28901.SAMN13057273.GCA_009230165_01317': 'no_eggnog'}

    Returns:
    - Dictionary of aggregated statistics for all batches {func: {mge_type: {cluster: COG_category}}}
    """
    aggregated = {}

    for batch in batches:
        for func_name, mge_dict in batch.items():
            aggregated[func_name] = {}
            for mge_type, attr_dict in mge_dict.items():
                if mge_type not in aggregated[func_name]:
                    aggregated[func_name][mge_type] = defaultdict(list)  # TODO: Generalise to function 
                # Update using the contents of the nested dictionary
                for cluster_id, value in attr_dict.items():
                    aggregated[func_name][mge_type][cluster_id] = value # TODO: replace with proper majority vote 
    # flatten the COG value per batch 
    return aggregated

def aggregate_counts(batch_counts):
    """
    Aggregate statistics across all batches.

    Parameters:
    - batch_counts: List of batch statistics (dictionaries with Counter values).

    Returns:
    - Dictionary of aggregated statistics for all batches {func: {mge_type: {cluster: count}}}
    """
    aggregated = {}

    for batch in batch_counts:
        for func_name, mge_counter in batch.items():
            aggregated[func_name] = {}
            for mge_type, nested_counter in mge_counter.items():
                if mge_type not in aggregated[func_name]:
                    aggregated[func_name][mge_type] = Counter()  # Initialize if not already present
                # Update using the contents of the nested dictionary
                for cluster_id, value in nested_counter.items():
                    if isinstance(value, str):
                        aggregated[func_name][mge_type][cluster_id] = value # Overwrite previous COG category with the latest batch 
                    else:
                        if cluster_id in aggregated[func_name][mge_type]:
                            aggregated[func_name][mge_type][cluster_id] += value
                        else:
                            aggregated[func_name][mge_type][cluster_id] = value
    return aggregated

