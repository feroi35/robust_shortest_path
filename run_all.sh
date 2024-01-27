method="dualized"
time_limit=600
verbose=0

GLOBAL_CSV="results/${method}_results_ter.csv"
repo_instances="data/processed"

# Supprimer le fichier CSV global s'il existe déjà
rm -f $GLOBAL_CSV

# En-tête du fichier CSV global
if method="heuristics"; then
    echo "instance,method,objective,unproved_lower_bound,time,nb_nodes,robust_constraint,static_objective,static_constraint,S,path" > $GLOBAL_CSV
else
    echo "instance,method,objective,lower_bound,time,nb_nodes,robust_constraint,static_objective,static_constraint,S,path" > $GLOBAL_CSV
fi

make release

count=1
for file in "$repo_instances"/*.gr
do
    echo "Running file: $file, $count/123"
    ./myprogram "$file" "$method" "$time_limit" "$verbose" >> "$GLOBAL_CSV"
    ((count++))
    # if [ "$count" -lt 5 ]; then
    #     echo "Running file: $file, $count/123"
    #     ./myprogram "$file" "$method" $verbose >> "$GLOBAL_CSV"
    #     ((count++))
    # else
    #     break
    # fi
done

# make clean
