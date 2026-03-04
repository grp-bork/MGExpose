# pylint: disable=C0301,E0401,W1510,W0621
# flake8: noqa

'''
mge_annotation.py
GCA_000012825.1.genomes
GCA_000012825.1.genomes.gff.gz
GCA_000012825.1.genomes.recombinase_hmmsearch.besthits.out
specI_v4_00061
txsscan_rules.txt
all_systems.tsv
GCA_000012825.1.genomes.emapper.annotations
mge_rules_ms.txt
--cluster_data GCA_000012825.1.genomes_mmseqcluster.tsv.gz
--output_dir specI_v4_00061/GCA_000012825.1.genomes/
--dump_intermediate_steps
'''
import os
import subprocess
import pytest

TEST_DATADIR = "../test_data/current/"
OUTPUT_DIR = "specI_v4_00061/GCA_000012825.1.genomes/"
GFF_FILENAME = "GCA_000012825.1.genomes.full_length_MGE_assignments.gff3"
TXT_FILENAME = "GCA_000012825.1.genomes.full_length_MGE_assignments.txt"
DEBUG_FILES = {
    # "GCA_000012825.1.genomes.assign_mge.step1.txt": "GCA_000012825.1.genomes.assign_mge.step1.txt",
    # "GCA_000012825.1.genomes.assign_mge.step2.txt": "GCA_000012825.1.genomes.assign_mge.step2.txt",
    # "GCA_000012825.1.genomes.assign_mge.step3.txt": "GCA_000012825.1.genomes.assign_mge.step3.txt",
    # "GCA_000012825.1.genomes.pan_genome_calls.txt": "GCA_000012825.1.genomes.pan_genome_calls.txt",
    # "GCA_000012825.1.genomes.pan_genome_islands.txt": "GCA_000012825.1.genomes.pan_genome_islands.txt",
}

TEST_OUT_GFF = os.path.join(TEST_DATADIR, "output", OUTPUT_DIR, GFF_FILENAME)
TEST_OUT_TXT = os.path.join(TEST_DATADIR, "output", OUTPUT_DIR, TXT_FILENAME)

INPUT_ARGS = [
    ("genome_id", "GCA_000012825.1.genomes"),
    ("prodigal_gff", os.path.join(TEST_DATADIR, "GCA_000012825.1.genomes.gff.gz")),
    ("recombinase_hits", os.path.join(TEST_DATADIR, "GCA_000012825.1.genomes.recombinase_hmmsearch.besthits.out")),
    ("mge_rules", os.path.join(TEST_DATADIR, "mge_rules_ms.txt")),
    ("--speci", "specI_v4_00061"),
    ("--txs_macsy_rules", os.path.join(TEST_DATADIR, "txsscan_rules.txt")),
    ("--txs_macsy_report", os.path.join(TEST_DATADIR, "all_systems.tsv")),
    ("--phage_eggnog_data", os.path.join(TEST_DATADIR, "GCA_000012825.1.genomes.emapper.annotations")),
    ("--cluster_data", os.path.join(TEST_DATADIR, "GCA_000012825.1.genomes_mmseqcluster.tsv.gz")),
    ("--output_dir", OUTPUT_DIR),
    ("--dump_intermediate_steps", ""),
]


@pytest.fixture(scope="module")
def run_mge_annotation(tmpdir_factory):
    """Fixture to run the mge_annotation.py script once per module and generate the output files."""
    tmpdir = tmpdir_factory.mktemp("mge_output")
    tmp_output_dir = os.path.join(tmpdir, OUTPUT_DIR)
    os.makedirs(tmp_output_dir, exist_ok=True)
    debug_dir = None

    # Prepare the command with input arguments
    command = ["python", "mge_annotation.py", "denovo"]
    for arg, val in INPUT_ARGS:
        if val:  # Append argument only if value is non-empty
            if arg == "--output_dir":
                command.extend(["--output_dir", tmp_output_dir])
            elif '--' in arg:
                command.extend([arg, val])
            else:
                command.append(val)  # Obligatory input
        elif arg == "--dump_intermediate_steps":
            command.append(arg)
            debug_dir = os.path.join(tmp_output_dir, "debug")
            os.makedirs(debug_dir, exist_ok=True)

        else:
            command.append(arg)

    # Execute the command
    print("Running command:", command)
    result = subprocess.run(command, capture_output=True, text=True)

    # Ensure the script ran successfully
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"
    return tmp_output_dir, debug_dir


def compare_output_files(generated_file_path, expected_file_path):
    """Helper function to compare a generated file with its expected output."""
    # Ensure the generated file exists
    assert os.path.exists(generated_file_path), f"Generated file not found at {generated_file_path}"

    # Read both the expected file and the generated file
    with open(expected_file_path, "rb") as f:
        expected_content = f.read()

    with open(generated_file_path, "rb") as f:
        generated_content = f.read()

    # Assert that the content of both files is identical
    assert expected_content == generated_content, f"The generated file {generated_file_path} does not match the expected output."


def test_gff_output(run_mge_annotation):
    """Test to compare the generated GFF file with the expected output."""
    tmp_output_dir, _ = run_mge_annotation
    print("Temporary output directory: ", tmp_output_dir)
    generated_gff_path = os.path.join(tmp_output_dir, GFF_FILENAME)
    compare_output_files(generated_gff_path, TEST_OUT_GFF)


def test_txt_output(run_mge_annotation):
    """Test to compare the generated TXT file with the expected output."""
    tmp_output_dir, _ = run_mge_annotation
    generated_txt_path = os.path.join(tmp_output_dir, TXT_FILENAME)
    compare_output_files(generated_txt_path, TEST_OUT_TXT)


# Individual tests for each debug file
@pytest.mark.parametrize("debug_filename", DEBUG_FILES.keys())
def test_debug_file_output(run_mge_annotation, debug_filename):
    """Test to compare each file in the debug directory with the expected output."""
    _, debug_dir = run_mge_annotation
    generated_file_path = os.path.join(debug_dir, debug_filename)
    expected_file_path = os.path.join(TEST_DATADIR, "output", OUTPUT_DIR, "debug", debug_filename)
    compare_output_files(generated_file_path, expected_file_path)
