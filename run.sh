GLOBAL_CSV="results/static_results.csv"
repo_instances="data/processed"

# Supprimer le fichier CSV global s'il existe déjà
rm -f $GLOBAL_CSV

# En-tête du fichier CSV global
echo "instance,objective,time" > $GLOBAL_CSV

make release

count=0
for file in "$repo_instances"/*.gr
do
    echo "Running file: $file"
    ./myprogram "$file" >> "$GLOBAL_CSV"
    # if [ "$count" -lt 2 ]; then
    #     echo "Running file: $file"
    #     ./myprogram "$file" >> "$GLOBAL_CSV"
    #     ((count++))
    # else
    #     break
    # fi
done

make clean
