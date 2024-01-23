# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"


method="heuristics"
verbose=0


make release

echo "Running file: $file"
./myprogram "$file" "$method" "$verbose"

# make clean
