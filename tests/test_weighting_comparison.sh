#!/bin/bash
# Test all weighting types on EIS_ZRY-3d-1.dta

echo "======================================================================="
echo "Srovnání weighting typů na problematických datech EIS_ZRY-3d-1.dta"
echo "======================================================================="
echo ""

for wtype in uniform sqrt proportional square; do
    echo "-----------------------------------------------------------------------"
    echo "Testing: $wtype"
    echo "-----------------------------------------------------------------------"

    output=$(./eis.py EIS_ZRY-3d-1.dta --voigt-chain --f-max 10 --voigt-n-per-decade 3 --voigt-prune-threshold 0.005 --weighting "$wtype" --no-show 2>&1)

    # Extract key info
    echo "$output" | grep -E "Vážení:" | head -1
    echo "$output" | grep -E "Initial guess summary:" -A 5 | grep "Voigt"
    echo "$output" | grep -E "Fit (error|chyba)" | tail -1

    # Look for the final line with error
    last_line=$(echo "$output" | tail -1)
    if [[ $last_line =~ ^[0-9]+ ]]; then
        error=$(echo "$last_line" | awk '{print $3}')
        echo "Final error: $error %"
    fi

    # Check for warnings
    if echo "$output" | grep -q "VYSOKÁ CHYBA"; then
        echo "⚠️  HIGH ERROR WARNING detected"
    fi
    if echo "$output" | grep -q "Vysoká korelace"; then
        echo "⚠️  HIGH CORRELATION WARNING detected"
    fi
    if echo "$output" | grep -q "vysoká uncertainty"; then
        echo "⚠️  HIGH UNCERTAINTY WARNING detected"
    fi

    echo ""
done

echo "======================================================================="
echo "Test complete"
echo "======================================================================="
