   user  system elapsed 
990.030 375.069 220.588 


base R, No opt, linked with openblas: 
   user.self sys.self elapsed user.child sys.child
t1   643.363   32.314 122.944          0         0
t2  2194.482 1562.861 537.891          0         0

base R, with opt, linked with openblas: 
   user.self sys.self elapsed user.child sys.child
t1   666.381   23.706 127.471          0         0
t2  2160.767 1574.761 514.817          0         0


base R, No opt, linked with MKL: 
base R, with opt, linked with MKL (static): 
   user.self sys.self elapsed user.child sys.child
t1   293.223    4.899 128.611          0         0
t2   547.201   18.472 110.949          0         0


base R, with opt, linked with MKL (dynamic): 
   user.self sys.self elapsed user.child sys.child
t1   287.014    6.017 127.732          0         0
t2   533.650   16.634 111.246          0         0


MKL R, No opt, MKL: 
MKL R, with opt, MKL: 


Mac OS with Accelerate
user.self sys.self elapsed user.child sys.child 
t1 1271.351 21.503 142.559 0 0 
t2 515.349 20.305 76.213 0 0 


require(ACTIONet)
download.file('http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5', 'pbmc_10k_v3.h5') 
ace = import.ace.from.10X.h5('pbmc_10k_v3.h5', prefilter = T, min_cells_per_feat = 0.01, min_feats_per_cell = 1000)
t1 = system.time( {ace = reduce.ace(ace, SVD_algorithm=0)} )
readr::write_rds(ace, "pbmc_10k_v3.rds")

#S_r = ACTIONet::colMaps(ace)[["ACTION"]]
#S = counts(ace)

require(ACTIONet)
ace = readr::read_rds("pbmc_10k_v3.rds")

t2 = system.time( {ace = run.ACTIONet(ace, min.cells.per.arch = 2, unification.resolution = 1)} )
#print(t2)
print(rbind(t1, t2))


data("curatedMarkers_human")
markers = curatedMarkers_human$Blood$PBMC$Ding2019$marker.genes
annot.out = annotate.cells.using.markers(ace, markers)
ace$celltypes = annot.out$Labels

plot.ACTIONet(ace, "celltypes", transparency.attr = ace$node_centrality)



library(jointprof)

# Profile now
Sys.setenv(PPROF_BINARY_PATH="~/go/bin/pprof")
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/usr/local/go/bin/", sep = ":"))
out_file <- "./ACTION_jointprof.out"
start_profiler(out_file)
#ACTION.out = run_ACTION(S_r, k_min = 2, k_max = 30, thread_no = 8, max_it = 50, min_delta = 1e-300)
reduction.out = reduce_kernel(S, reduced_dim = 50, iter = 5, seed = 0, reduction_algorithm = 1, SVD_algorithm = 1)
profile_data <- stop_profiler()

pprof_file <- tempfile("jointprof", fileext = ".pb.gz")
profile::write_pprof(profile_data, pprof_file)

## knit
dir.create(knit_dir("jointprof_fig"), recursive = TRUE)
svg_file <- knit_dir("jointprof_fig/minimal.svg")

system2(
  find_pprof(),
  c(
    "-svg",
    "-nodefraction 0.01",
    "-output",
    shQuote(svg_file),
    shQuote(pprof_file)
  )
)

png_file <- knit_dir("jointprof_fig/minimal.png")
rsvg::rsvg_png(svg_file, png_file, width = 680)









target_file <- "Rprof.out"
start_profiler(target_file)
reduction.out = reduce_kernel(S, reduced_dim = 50, iter = 5, seed = 0, reduction_algorithm = 1, SVD_algorithm = 1)
stop_profiler()









  #########################################
  #########################################
  #########################################
  path = "Rprof.out"
  numfiles = 100L
  bufsize = 10000L
  method = "run_ACTION"
  
  find_pprof()

  pprof_path <- sprintf("%s_jointprof.prof", method)
  rprof_path <- sprintf("%s_jointprof.out", method)

  prof_data <- jointprof:::init_profiler_impl()
  utils::Rprof(
	filename = rprof_path,
	line.profiling = TRUE,
	gc.profiling = TRUE,
	numfiles = numfiles,
	bufsize = bufsize
  )
  jointprof:::start_profiler_impl(prof_data, pprof_path)
  jointprof:::.my_env$prof_data <- prof_data
  jointprof:::.my_env$path <- path
  jointprof:::.my_env$pprof_path <- pprof_path
  jointprof:::.my_env$rprof_path <- rprof_path
  #########################################
  #########################################
  #########################################


reduction.out = reduce_kernel(S, reduced_dim = 50, iter = 5, seed = 0, reduction_algorithm = 1, SVD_algorithm = 1)


  #########################################
  #########################################
  #########################################
  utils::Rprof(NULL)
  jointprof:::stop_profiler_impl(jointprof:::.my_env$prof_data)

  ds <- jointprof:::combine_profiles(.jointprof:::my_env$pprof_path, jointprof:::.my_env$rprof_path)
  profile::write_rprof(ds, .my_env$path)
  #########################################
  #########################################
  #########################################








