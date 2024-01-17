# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"

# debuguer avec 1200-COL 1500-NY 1400-COL tout les 2000: probablement la time limit permet pas de resoudre la racine

# method can be:
#   - static (for static problem)
#   - dualized (for dualized formulation)
method="dualized"

# verbose
#   - 0 (no verbose)
#   - 1 (minimum verbose)
#   - 2 (maximum verbose)
verbose=1

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$verbose"

# make clean
