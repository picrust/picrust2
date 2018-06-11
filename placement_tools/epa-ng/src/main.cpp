#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include <CLI/CLI.hpp>

#include "net/mpihead.hpp"
#include "util/logging.hpp"
#include "util/Options.hpp"
#include "util/stringify.hpp"
#include "util/parse_model.hpp"
#include "util/split.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/Binary.hpp"
#include "io/file_io.hpp"
#include "io/msa_reader.hpp"
#include "tree/Tree.hpp"
#include "core/raxml/Model.hpp"
#include "core/place.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"

#ifndef EPA_VERSION
#define EPA_VERSION "UNKNOWN"
#endif

static void ensure_dir_has_slash(std::string& dir)
{
  if (dir.length() > 0 && dir.back() != '/') {
    dir += "/";
  }
}

inline bool is_file (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

void exit_epa(int ret=EXIT_SUCCESS)
{
  MPI_FINALIZE();
  std::exit(ret);
}

int main(int argc, char** argv)
{
  auto start_all = std::chrono::high_resolution_clock::now();
  genesis::utils::Logging::log_to_stdout();

#ifdef __MPI
  MPI_INIT(&argc, &argv);
  int local_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  if (local_rank != 0) {
    genesis::utils::Logging::log_to_stdout(false);  
  }
#endif
  genesis::utils::Logging::max_level(genesis::utils::Logging::kInfo);

  std::string invocation("");
  std::string model_desc("GTR+G");
  Options options;

  for (int i = 0; i < argc; ++i) {
    invocation += argv[i];
    invocation += " ";
  }

  std::string query_file;
  std::string work_dir("./");
  std::string tree_file;
  std::string reference_file;
  std::string binary_file;
  std::string bfast_conv_file;
  std::vector<std::string> split_files;

  std::string banner;

  raxml::Model model;

  bool display_version  = false;
  bool verbosity        = false;
  bool heuristics_off   = not options.prescoring;
  bool raxml_blo        = not options.sliding_blo;
  bool no_pre_mask      = not options.premasking;

  const bool empty = argc == 1;

  CLI::App app{"epa-ng - Massively-Parallel Evolutionary Placement Algorithm"};

  //  ============== GENERAL OPTIONS ==============

  app.add_flag("-v,--version", display_version, "Display version.");
  app.add_flag("--verbose", verbosity, "Display debug output.");

  //  ============== CONVERT OPTIONS ==============

  app.add_option( "-c,--bfast",
                  bfast_conv_file,
                  "Convert the given fasta file to bfast format."
                )->group("Convert")->check(CLI::ExistingFile);
  app.add_flag( "-B,--dump-binary",
                  options.dump_binary_mode,
                  "Binary Dump mode: write ref. tree in binary format then exit. NOTE: not compatible with premasking!"
                )->group("Convert");
  app.add_option( "--split",
                  split_files,
                  "Takes a reference MSA (phylip) and combined ref +"
                  " query MSA(s) (phylip) and outputs one pure query file (fasta). "
                  "Usage: epa-ng --split ref_alignment query_alignments+"
                )->group("Convert")->check(CLI::ExistingFile);

  //  ============== INPUT OPTIONS ==============
  auto tree_file_opt =
  app.add_option( "-t,--tree",
                  tree_file,
                  "Path to Reference Tree file."
                )->group("Input")->check(CLI::ExistingFile);
  auto reference_file_opt =
  app.add_option( "-s,--ref-msa",
                  reference_file,
                  "Path to Reference MSA file."
                )->group("Input")->check(CLI::ExistingFile);
  auto binary_file_opt =
  app.add_option( "-b,--binary",
                  binary_file,
                  "Path to binary reference file, as created using --dump-binary."
                )->group("Input")->check(CLI::ExistingFile);

  binary_file_opt->excludes(tree_file_opt)->excludes(reference_file_opt);
  tree_file_opt->excludes(binary_file_opt);  
  reference_file_opt->excludes(binary_file_opt);  

  app.add_option( "-q,--query",
                  query_file,
                  "Path to Query MSA file."
                )->group("Input")->check(CLI::ExistingFile);
  app.add_option( "-m,--model",
                  model_desc,
                  "Description string of the model to be used, or a RAxML_info file."
                  " --model STRING | FILE "
                  "See: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model",
                  true
                )->group("Input");

  //  ============== OUTPUT OPTIONS ==============

  app.add_option("-w,--outdir", work_dir, "Path to output directory.", true
                )->group("Output")->check(CLI::ExistingDirectory);

  app.add_option("--tmp", options.tmp_dir, "Path to temporary directory. If set, MPI-Rank-local"
                  " files will be stored here instead. Useful for node-local SSDs!"
              )->group("")->check(CLI::ExistingDirectory)
  // ->set_custom_option("DIR", 2)
  ;
  auto filter_acc_lwr =
  app.add_option( "--filter-acc-lwr",
                  options.support_threshold,
                  "Accumulated likelihood weight after which further placements are discarded.",
                  options.acc_threshold
                )->group("Output")->check(CLI::Range(0.0,1.0));
  auto filter_min_lwr =
  app.add_option( "--filter-min-lwr",
                  options.support_threshold,
                  "Minimum likelihood weight below which a placement is discarded.",
                  not options.acc_threshold
                )->group("Output")->check(CLI::Range(0.0,1.0));
  filter_acc_lwr->excludes(filter_min_lwr);
  filter_min_lwr->excludes(filter_acc_lwr);

  auto filter_min =
  app.add_option( "--filter-min",
                  options.filter_min,
                  "Minimum number of placements per sequence to include in final output.",
                  true
                )->group("Output");
  auto filter_max =
  app.add_option( "--filter-max",
                  options.filter_max,
                  "Maximum number of placements per sequence to include in final output.",
                  true
                )->group("Output");

  //  ============== COMPUTE OPTIONS ==============

  auto dyn_heur =
  app.add_option( "-g,--dyn-heur",
                  options.prescoring_threshold,
                  "Two-phase heuristic, determination of candidate edges using accumulative threshold. "
                  "Enabled by default! See --no-heur for disabling it",
                  true
                )->group("Compute")->check(CLI::Range(0.0,1.0));
  auto fix_heur =
  app.add_option( "-G,--fix-heur",
                  options.prescoring_threshold,
                  "Two-phase heuristic, determination of candidate edges by specified percentage of total edges.",
                  options.prescoring_by_percentage
                )->group("Compute")->check(CLI::Range(0.0,1.0));
  auto baseball_heur =
  app.add_flag( "--baseball-heur",
                  options.baseball,
                  "Baseball heuristic as known from pplacer. strike_box=3,max_strikes=6,max_pitches=40."
                )->group("Compute");
  auto no_heur =
  app.add_flag( "--no-heur",
                  heuristics_off,
                  "Disables heuristic preplacement completely. Overrides all other heuristic flags."
                )->group("Compute");
  dyn_heur->excludes(fix_heur)->excludes(baseball_heur)->excludes(no_heur);
  fix_heur->excludes(dyn_heur)->excludes(baseball_heur)->excludes(no_heur);
  baseball_heur->excludes(dyn_heur)->excludes(fix_heur)->excludes(no_heur);
  no_heur->excludes(dyn_heur)->excludes(fix_heur)->excludes(baseball_heur);

  auto chunk_size =
  app.add_option( "--chunk-size",
                  options.chunk_size,
                  "Number of query sequences to be read in at a time. May influence performance.",
                  true
                )->group("Compute");
  app.add_flag( "--raxml-blo",
                  raxml_blo,
                  "Employ old style of branch length optimization during thorough insertion as opposed to sliding approach. "
                  "WARNING: may significantly slow down computation."
                )->group("Compute");
  app.add_flag( "--no-pre-mask",
                  no_pre_mask,
                  "Do NOT pre-mask sequences. Enables repeats unless --no-repeats is also specified."
                )->group("Compute");

  #ifdef __OMP
  auto threads =
  app.add_option( "-T,--threads",
                  options.num_threads,
                  "Number of threads to use. If 0 is passed as argument,"
                  "program will run with the maximum number of threads available.",
                  true
                )->group("Compute");
  #endif

  try {
    app.parse(argc, argv);
    if (empty) {
      throw CLI::CallForHelp();
    }
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  if (display_version) {
    std::cout << "EPA-ng v" << EPA_VERSION << std::endl;
    exit_epa();
  }

  ensure_dir_has_slash(work_dir);
  if ( not options.tmp_dir.empty() ) {
    ensure_dir_has_slash(options.tmp_dir);
    LOG_INFO << "Selected: Temporary dir: " << options.tmp_dir;
  }

  // no log file for conversion functions
  if (not bfast_conv_file.empty()) {
    LOG_INFO << "Converting given FASTA file to BFAST format...";
    auto resultfile = Binary_Fasta::fasta_to_bfast(bfast_conv_file, work_dir);
    LOG_INFO << "Resulting bfast file was written to: " << resultfile;
    exit_epa();
  }

  if (split_files.size()) {
    if (split_files.size() < 2) {
      LOG_ERR << "Incorrect number of inputs! Usage: epa-ng --split ref_alignment query_alignments+";
      exit_epa();
    }
    auto ref_msa = split_files[0];
    split_files.erase(split_files.begin());
    LOG_INFO << "Splitting files based on reference: " << ref_msa;
    split(ref_msa, split_files, work_dir);
    exit_epa();
  }

  #ifdef __MPI
  genesis::utils::Logging::log_to_file(work_dir + std::to_string(local_rank) + ".epa_info.log");
  #else
  genesis::utils::Logging::log_to_file(work_dir + "epa_info.log");
  #endif

  LOG_INFO << "Selected: Output dir: " << work_dir;

  if (verbosity) {
    LOG_INFO << "Selected: verbose (debug) output";
    genesis::utils::Logging::max_level(genesis::utils::Logging::kDebug);
  }
  
  if (not query_file.empty()) {
    LOG_INFO << "Selected: Query file: " << query_file;
  }

  if (not tree_file.empty()) {
    LOG_INFO << "Selected: Tree file: " << tree_file;
  }

  if (not reference_file.empty()) {
    LOG_INFO << "Selected: Reference MSA: " << reference_file;
  }

  if (not binary_file.empty()) {
    options.load_binary_mode = true;
    LOG_INFO << "Selected: Binary CLV store: " << binary_file;
  }

  if (*filter_acc_lwr)
  {
    options.acc_threshold = true;
    LOG_INFO << "Selected: Filtering by accumulated threshold: " << options.support_threshold;
  }

  if (*filter_min_lwr) {
    options.acc_threshold = false;
    LOG_INFO << "Selected: Filtering by minimum threshold: " << options.support_threshold;
  }

  if (*filter_min) {
    LOG_INFO << "Selected: Minimum number of placements per query: " << options.filter_min;
  }

  if (*filter_max) {
    LOG_INFO << "Selected: Maximum number of placements per query: " << options.filter_max;
  }

  if (options.filter_min > options.filter_max) {
    throw std::runtime_error{"filter-min must not exceed filter-max!"};
  }

  if (*fix_heur) {
    options.prescoring = options.prescoring_by_percentage = true;
    LOG_INFO << "Selected: Prescoring by percentage of branches: " << options.prescoring_threshold;
  }

  if (*dyn_heur) {
    options.prescoring = true;
    options.prescoring_by_percentage = false;
    LOG_INFO << "Selected: Prescoring by accumulated LWR threshold: " << options.prescoring_threshold;
  }

  if (*baseball_heur) {
    LOG_INFO << "Selected: Prescoring using the baseball heuristic";
  }

  if (raxml_blo) {
    options.sliding_blo = false;
    LOG_INFO << "Selected: On query insertion, optimize branch lengths the way RAxML-EPA did it";
  }

  if (no_pre_mask) {
    options.premasking = false;
    options.repeats = true;
    LOG_INFO << "Selected: Disabling pre-masking. (repeats enabled!)";
  }

  // if (cli.count("no-repeats")) {
  //   options.repeats = false;
  //   LOG_INFO << "Selected: Using the non-repeats version of libpll/modules";
  // }

  if (*no_heur) {
    options.prescoring = false;
    LOG_INFO << "Selected: Disabling the prescoring heuristics.";
  }

  if (options.dump_binary_mode) {
    LOG_INFO << "Selected: Build reference tree and write it out as a binary CLV store (for MPI)";
    LOG_INFO << "\tWARNING: this mode means that no placement will take place in this run";
  }

  if (is_file(model_desc)) {
    LOG_INFO << "Selected: Specified model file: " << model_desc;
    model_desc = parse_model(model_desc);
  } else {
    LOG_INFO << "Selected: Specified model: " << model_desc;
  }

  model = raxml::Model(model_desc);

  LOG_INFO << model;

  if (*chunk_size) {
    LOG_INFO << "Selected: Reading queries in chunks of: " << options.chunk_size;
  }
  #ifdef __OMP
  if (*threads) {
    LOG_INFO << "Selected: Using threads: " << options.num_threads;
  }
  #endif

  //================================================================
  //============    EPA    =========================================
  //================================================================

  banner += "    ______ ____   ___           _   __ ______\n";
  banner += "   / ____// __ \\ /   |         / | / // ____/\n";
  banner += "  / __/  / /_/ // /| | ______ /  |/ // / __  \n";
  banner += " / /___ / ____// ___ |/_____// /|  // /_/ /  \n";
  banner += "/_____//_/    /_/  |_|      /_/ |_/ \\____/ (v";
  banner += EPA_VERSION;
  banner += ")\n \n";

  LOG_INFO << banner << std::endl;

  LOG_DBG << "Peeking into MSA files and generating masks";

  MSA_Info ref_info;
  if (not reference_file.empty()) {
    ref_info = make_msa_info(reference_file);
    LOG_DBG << "Reference File:\n" << ref_info;
  }

  MSA_Info qry_info;
  if (not query_file.empty()) {
    qry_info = make_msa_info(query_file);
    LOG_DBG << "Query File:\n" << qry_info;
  }

  MSA_Info::or_mask(ref_info, qry_info);

  // msa_info.reset_gaps <- --no-pre-mask

  MSA ref_msa;
  if (reference_file.size()) {
    ref_msa = build_MSA_from_file(reference_file, ref_info, options.premasking);
    LOG_DBG << "Reference File size: " << ref_msa.size();
    LOG_DBG << "Reference File width: " << ref_msa.num_sites();
    if (ref_msa.size() == 0 or ref_msa.num_sites() == 0 ) {
      throw std::runtime_error{"Something went wrong reading the reference file."};
    }
  }

  // build the Tree
  Tree tree;
  if (options.load_binary_mode) {
    LOG_INFO << "Loading from binary";
    tree = Tree(binary_file, model, options);
  } else {
    // build the full tree with all possible clv's
    tree = Tree(tree_file, ref_msa, model, options);
  }

  if (not options.dump_binary_mode) {
    if (query_file.empty()) {
      throw std::runtime_error{"Must supply query file! Combined MSA files not currently supported, please split them and specify using -s and -q."};
    }
  } else {
    // dump to binary if specified
    LOG_INFO << "Writing to binary";
    std::string dump_file(work_dir + "epa_binary_file");
    dump_to_binary(tree, dump_file);
    exit_epa();
  }

  // start the placement process and write to file
  auto start_place = std::chrono::high_resolution_clock::now();
  simple_mpi(tree, query_file, qry_info, work_dir, options, invocation);
  auto end_place = std::chrono::high_resolution_clock::now();
  auto placetime = std::chrono::duration_cast<std::chrono::seconds>(end_place - start_place).count();

  LOG_INFO << "Time spent placing: " << placetime << "s";

  MPI_FINALIZE();

  auto end_all = std::chrono::high_resolution_clock::now();
  auto alltime = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all).count();

  LOG_INFO << "Elapsed Time: " << alltime << "s";

	return EXIT_SUCCESS;
}
