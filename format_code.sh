# clang-format -style=google -dump-config > .clang-format
# find . \( -path './inst/include/ACTIONet/leiden/*' -o -path './src/ACTIONet/network_tools/leiden/*' -o -path './inst/include/ACTIONet/boost/*' \) -prune -o -type f -regex '.*\.\(cc\|c\|h\)' -exec clang-format -style=file -i {} \;
# Rscript -e "formatR::tidy_dir('`pwd`/R', recursive = T)"
