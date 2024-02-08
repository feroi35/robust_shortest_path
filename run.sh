

file="data/processed/20_USA-road-d.COL.gr"
method="plans_coupants"
time_limit=50
verbose=2

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$time_limit" "$verbose"

make clean
