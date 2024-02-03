file="data/processed/1000_USA-road-d.COL.gr"
# file="data/processed/20_USA-road-d.BAY.gr"
# file="data/retest_dualized/1100_USA-road-d.BAY.gr"

method="dualized"
time_limit=500

verbose=2

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

# make clean
