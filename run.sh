# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"

make release

echo "Running file: $file"
./myprogram "$file"

# make clean
