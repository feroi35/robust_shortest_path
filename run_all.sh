method="dualized"
verbose="0"

# faudrait essayer que en changeant la methode ca change le fichier de sortie
GLOBAL_CSV="results/dualized_results.csv"
repo_instances="data/processed"

# Supprimer le fichier CSV global s'il existe déjà
rm -f $GLOBAL_CSV

# En-tête du fichier CSV global
echo "instance,objective,time,nb_nodes,lower_bound" > $GLOBAL_CSV

make release

count=0
for file in "$repo_instances"/*.gr
do
    echo "Running file: $file"
    ./myprogram "$file" "$method" "$verbose" >> "$GLOBAL_CSV"
    # if [ "$count" -lt 2 ]; then
    #     echo "Running file: $file"
    #     ./myprogram "$file" >> "$GLOBAL_CSV"
    #     ((count++))
    # else
    #     break
    # fi
done

make clean
