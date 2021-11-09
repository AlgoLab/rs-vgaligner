#!/bin/bash 

# Install vgaligner (just for safety)
cargo install --path ../..

for d in */; do
	echo "\n"
	echo "========Executing cleaning pipeline for $d======="

	# Move to d
	cd "$d"

	# Remove txt
	find . -name "*.txt" -type f -delete

	# Remove gaf
	find . -name "*.gaf" -type f -delete

	# Remove odgi
	find . -name "*.odgi" -type f -delete

	# Move back to main folder
	cd ..

	echo "========End pipeline for $d======="
done
