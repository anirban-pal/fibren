#!/bin/bash

echo "ref"
cat $1 | grep "rel_error: 1.00e-15" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
echo "fibren_adaptive"
cat $1 | grep "rel_error: 1.00e-02" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-03" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-04" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-05" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-06" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-07" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-08" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-09" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "rel_error: 1.00e-10" | awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
echo "fibren_discrete"
cat $1 | grep "n: 1 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 3 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 10 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 31 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 100 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 316 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 1000 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 3162 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 10000 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
cat $1 | grep "n: 31622 " |  awk '{sum+=$8; sumsq+=$8*$8} END {print $0, sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}' | awk '{print $2" "$4" "$10" "$11}'
