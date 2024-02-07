
# file="data/processed/20_USA-road-d.COL.gr"
file="data/processed/1000_USA-road-d.COL.gr"
# file="data/processed/1200_USA-road-d.BAY.gr"
# file="data/processed/1300_USA-road-d.BAY.gr"

# method="static"
# method="heuristic"
# method="branch_and_cut"
# method="dualized"
method="branch_and_cut"
time_limit=10

verbose=0

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

# make clean
