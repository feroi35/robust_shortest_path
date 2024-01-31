# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"

method="static"
time_limit=600
verbose=2

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

# make clean
