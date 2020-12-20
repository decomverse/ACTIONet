clang-format -style=google -dump-config > .clang-format
find . \( -path './include/ACTIONet/leiden/*' -o -path './src/ACTIONet/network_tools/leiden/*' -o -path './include/ACTIONet/boost/*' \) -prune -o -type f -regex '.*\.\(cc\|c\|h\)' -exec clang-format -style=file -i {} \;
black python_interface/
