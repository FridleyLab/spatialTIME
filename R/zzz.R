.onAttach = function(...){
  packageStartupMessage("spatialTIME version:\n", crayon::white(packageVersion("spatialTIME")),
                        "\nIf using for publication, please cite our manuscript:\n",
                        crayon::white("https://doi.org/10.1093/bioinformatics/btab757"))
}
