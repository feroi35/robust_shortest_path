
file="data/processed/20_USA-road-d.COL.gr"
method="branch_and_cut"
time_limit=60
verbose=2

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

make clean
