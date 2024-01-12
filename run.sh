# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"

# method can be:
#   - static (for static problem)
#   - dualized (for dualized formulation)
method="dualized"

make release

echo "Running file: $file"
./myprogram "$file" "$method"

# make clean
