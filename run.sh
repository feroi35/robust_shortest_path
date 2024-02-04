
# file="data/processed/20_USA-road-d.COL.gr"
file="data/processed/2000_USA-road-d.BAY.gr"
# file="data/retest_dualized/1100_USA-road-d.BAY.gr"
# file="data/test/toy_symetries_removal.gr"

method="branch_and_cut"
time_limit=300


verbose=0

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

# make clean
