rm ../docs/CHANGELOG.rst
echo "*********" > ../docs/CHANGELOG.rst
echo "Changelog">> ../docs/CHANGELOG.rst
echo "*********" >> ../docs/CHANGELOG.rst
releases=("v2.0.3b" "v2.0.4" "v2.0.5" "v2.0.6b" "v2.0.7" "v2.0.8b" "v2.0.9" "v2.0.10" "v2.0.12" "v2.0.15" "v2.0.18" "v2.1.2" "v2.1.5" "v2.1.7")
num_elements=${#releases[@]}
for i in `seq 1 $((num_elements-1))`
do
    next_release=$i
    previous_release=$((i-1))

    python release_notes.py ${releases[$previous_release]} ${releases[$next_release]} >> ../docs/CHANGELOG.rst
    echo "">> ../docs/CHANGELOG.rst
done
