count=$1
x=0
echo "STATS:"
while [ $x -lt $count ]
    do
        echo "TEST $x STATS:"
        cat "test.sim.$x.stats"
        x=$(( $x + 1 ))
        echo
    done
echo "SUMMARY STATS:"
cat "summary.stats"
