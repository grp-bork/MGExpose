# pylint: disable=C0301,C0303
"""module docstring"""
import os
import subprocess


def run_emapper(args):
    """ docstring """

    cmd = f"emapper.py " \
            f"-i {args.proteins_fasta} " \
            f"--cpu {args.threads} " \
            f"--data_dir {args.eggnog_db} " \
            f"--dmnd_db {os.path.join(args.eggnog_db, 'eggnog_proteins.dmnd')} " \
            f"--output {os.path.join(args.output_dir, args.genome_id)} " \
            "-m diamond --dmnd_algo 0"

    with subprocess.Popen(
        cmd, shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ) as process:
        out, err = process.communicate()

        if process.returncode != 0:
            raise ValueError(
                f"Error: eggnog-mapper annotation failed with error code {process.returncode}:\n\n" 
                f"{out.decode()}\n"
                f"{err.decode()}\n"
            )

    print(out.decode().strip())

# def generate_args(emapper_input, output_prefix, output_dir, data_dir, scratch_dir, override, resume, n_cpus=1,):
#     return Namespace(
#         cpu = n_cpus,
#         cpus_per_worker = n_cpus,  # int(args.cpu / total_workers)
#         data_dir = data_dir,
#         dmnd_db = f"{data_dir}/eggnog_proteins.dmnd",
#         input = emapper_input,
#         output_dir = output_dir,
#         output = output_prefix,
#         override = override,
#         resume = resume,
#         scratch_dir = scratch_dir,
#         allow_overlaps = "none",
#         annotate_hits_table = None,
#         cache_file = None,
#         dbmem = False,
#         decorate_gff = "no",
#         decorate_gff_ID_field = "ID",
#         dmnd_algo = "0",
#         dmnd_block_size = None,
#         dmnd_evalue = 0.001,
#         dmnd_frameshift = None,
#         dmnd_ignore_warnings = True,
#         dmnd_index_chunks = None,
#         dmnd_iterate = "yes",
#         dmnd_score = None,
#         end_port = 53200,
#         evalue = 0.001,
#         excel = False,
#         excluded_taxa = None,
#         gapextend = None,
#         gapopen = None,
#         genepred = "search",
#         itype = "proteins",
#         go_evidence = "non-electronic",
#         go_excluded = {"ND", "IEA"},
#         matrix = None,
#         md5 = False,
#         mode = "diamond",
#         no_annot = False,
#         no_file_comments = False,
#         num_servers = 1,
#         num_workers = 1,
#         outfmt_short = False,
#         overlap_tol = 0.0,
#         pfam_realign = "none",
#         pident = None,
#         port = 51700,
#         query_cover = None,
#         report_orthologs = False,
#         score = None,
#         seed_ortholog_evalue = 0.001,
#         seed_ortholog_score = None,
#         sensmode = "sensitive",
#         subject_cover = None,
#         target_orthologs = "all",
#         target_taxa = None,
#         tax_scope_mode = "inner_narrowest",
#         tax_scope_ids = parse_tax_scope("auto"),
#         temp_dir = os.getcwd(),
#         timeout_load_server = 10,
#         trans_table = None,
#         translate = False,
#     )


# def run_emapper(emapper_input, output_prefix, data_dir, output_dir=".", scratch_dir=None, override=False, resume=False, n_cpus=1,):

#     # output_prefix = "emapper_out"
#     # output_dir = "."
#     # scratch_dir = None
#     # resume = False
#     # override = False
#     # emapper_input = sys.argv[1]
#     # emapper_diamond_db = "/scratch/schudoma/databases/eggnog_db/eggnog_proteins.dmnd"
#     # data_dir = "/scratch/schudoma/databases/eggnog_db"
#     # n_cpus = 8
#     try:
#         from eggnogmapper.annotation.tax_scopes.tax_scopes import parse_tax_scope
#         from eggnogmapper.common import get_citation, set_data_path
#         from eggnogmapper.emapper import Emapper
#         from eggnogmapper.emapperException import EmapperException
#     except ImportError:
#         print("Cannot run functional_annotation: eggnog-mapper is not installed.")
#         sys.exit(1)

#     set_data_path(data_dir)

#     try:
#         start_time = time.time()        

#         args = generate_args(
#             emapper_input, output_prefix, output_dir, data_dir, scratch_dir, override, resume, n_cpus=n_cpus,
#         )

#         emapper = Emapper(
#             args.itype, args.genepred, args.mode, (not args.no_annot),
#             args.excel, args.report_orthologs, args.decorate_gff,
#             args.output, args.output_dir, args.scratch_dir,
#             args.resume, args.override
#         )

#         # n, elapsed_time = emapper.run(args, args.input, args.annotate_hits_table, args.cache_file)
#         n, elapsed_time = emapper.run(args, args.input, args.annotate_hits_table, args.cache_file)
#         elapsed_time = time.time() - start_time

#         addons = [args.mode, args.genepred]
#         print(get_citation(addons))
#         print(f'Total hits processed: {n}')
#         print(f'Total time: {elapsed_time:.0f} secs')

#     except EmapperException as ee:
#         print(ee)
#         sys.exit(1)
#     except Exception as e:
#         traceback.print_exc()
#         # print(e)
#         sys.exit(1)
#     # else:
#     #     print("FINISHED")
#     #     sys.exit(0)

# #   emapper.py -i emapper.faa --data_dir ${eggnog_db} --output ${speci}/${genome_id}/${genome_id} -m diamond --cpu $task.cpus --dmnd_algo 0


#     # args.itype, args.genepred, args.mode, (not args.no_annot),
#     #                      args.excel, args.report_orthologs, args.decorate_gff,
#     #                      args.output, args.output_dir, args.scratch_dir,
#     #                      args.resume, args.override

#     #  def __init__(self, itype, genepred, mode, annot, excel, report_orthologs, decorate_gff, prefix, output_dir, scratch_dir, resume, override)
