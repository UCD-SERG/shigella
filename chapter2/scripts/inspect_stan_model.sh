# ==========================================================================
# inspect_stan_model.sh — Look at model_2.stan critical parts
#
# Goal: Verify the Stan model's definition of Omega_B matches what we
# expect (biomarker correlation, not parameter correlation).
#
# Run:
#   cd ~/chapter2
#   bash scripts/inspect_stan_model.sh 2>&1 | tee logs/inspect_stan.log
# ==========================================================================

set -euo pipefail

cd ~/chapter2

STAN=inst/stan/model_2.stan

if [ ! -f "$STAN" ]; then
    echo "ERROR: $STAN not found"
    exit 1
fi

echo "========================================================"
echo " INSPECT: model_2.stan critical parts"
echo "========================================================"
echo ""

# ----------------------------------------------------------------------
# 1. Parameters block — confirm L_Omega_B is K x K
# ----------------------------------------------------------------------
echo "=== Parameters block (where L_Omega_B is declared) ==="
sed -n '/^parameters {/,/^}/p' "$STAN"
echo ""

# ----------------------------------------------------------------------
# 2. Where is Omega_B computed in generated quantities?
# ----------------------------------------------------------------------
echo "=== generated quantities (where Omega_B is computed) ==="
sed -n '/^generated quantities {/,/^}/p' "$STAN" | head -40
echo ""

# ----------------------------------------------------------------------
# 3. Look for any obvious K vs P swap
# ----------------------------------------------------------------------
echo "=== Lines mentioning L_Omega_B, L_Omega_P, Omega_B, Omega_P ==="
grep -nE "L_Omega_B|L_Omega_P|Omega_B|Omega_P" "$STAN" || true
echo ""

# ----------------------------------------------------------------------
# 4. Kronecker function — order matters!
# ----------------------------------------------------------------------
echo "=== kron_chol function (Kronecker product) ==="
sed -n '/matrix kron_chol/,/^  }$/p' "$STAN"
echo ""

# ----------------------------------------------------------------------
# 5. How is Theta built from Z and the Kronecker factor?
# ----------------------------------------------------------------------
echo "=== Transformed parameters: Theta construction ==="
sed -n '/^transformed parameters {/,/^}/p' "$STAN"
echo ""

# ----------------------------------------------------------------------
# 6. K and P dimensions in data block
# ----------------------------------------------------------------------
echo "=== Data block ==="
sed -n '/^data {/,/^}/p' "$STAN"
echo ""

echo "========================================================"
echo " WHAT TO CHECK"
echo "========================================================"
echo ""
echo "  1. L_Omega_B should be cholesky_factor_corr[K]  (K = biomarkers)"
echo "  2. L_Omega_P should be cholesky_factor_corr[P]  (P = 5 kinetic params)"
echo "  3. kron_chol(L_Sigma_B, L_Sigma_P) — argument order matters"
echo "  4. vec(Theta) ordering: should match Kronecker output"
echo ""
echo " If L_Omega_B is on P (parameters) instead of K (biomarkers),"
echo " or if kron_chol arguments are swapped, then Omega_B[1,2]"
echo " is actually measuring something different than intended."
echo ""
