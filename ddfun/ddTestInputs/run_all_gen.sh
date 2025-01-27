# run all the generators

for file in *.py; do
   echo "Running $file"
   python $file
done