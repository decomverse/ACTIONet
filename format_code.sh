find -regex '.*/.*\.\(c\|cc\|cpp\|h\)$' | xargs clang-format -i
Rscript -e "formatR::tidy_dir('`pwd`/R', recursive = T)"
