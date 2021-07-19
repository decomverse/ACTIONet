clang-format -style=google -dump-config > .clang-format
find . \( -path './include/' -o -path './src/' \) -prune -o -type f -regex '.*\.\(cpp\|cc\|c\|h\|hpp\)' -exec clang-format -style=file -i {} \;
