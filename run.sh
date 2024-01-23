# file="data/processed/1200_USA-road-d.COL.gr"
file="data/processed/20_USA-road-d.BAY.gr"

<<<<<<< HEAD
method="heuristics"
verbose=0
=======
method="dualized"
verbose=2
>>>>>>> 51ca4d073226cdea3012f57bee78b355be008506

make release

echo "Running file: $file"
./myprogram "$file" "$method" "$verbose"

# make clean
