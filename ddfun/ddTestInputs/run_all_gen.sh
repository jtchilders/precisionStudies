# run all the generators

NTESTS=1000
for file in *.py; do
   echo "Running $file"
   python $file $NTESTS
done