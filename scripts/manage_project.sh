#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# --- Configuration ---
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$PROJECT_DIR/data"
REPORTS_DIR="$PROJECT_DIR/reports/ncar_eol_dee0099881"
DRAFTS_DIR="$DATA_DIR/drafts/figures"
SCRIPT_DIR="$PROJECT_DIR/scripts"

echo "=========================================================="
echo "   CASES-99 Pipeline Management & Environment Utility    "
echo "=========================================================="

# 1. Verify Julia Installation
if ! command -v julia &> /dev/null; then
    echo "❌ Error: Julia is not installed or not in your PATH."
    echo "Please install Julia before running this script."
    exit 1
fi
echo "✓ Julia executable detected: $(julia --version)"

# 2. Synchronize Project Environment & Dependencies
echo "🔄 Synchronizing Julia project environment status..."
julia --project="$PROJECT_DIR" -e '
    using Pkg
    using UUIDs  # Added to resolve the UUID function name error
    println("Checking registry status and instantiating environment...")
    Pkg.activate(".")
    Pkg.instantiate()

    # Explicitly ensure LaTeXStrings is added for the updated plotting script
    if !haskey(Pkg.dependencies(), UUID("b964fa9f-0449-5b57-a5c2-d326f3787e93"))
        println("Adding missing dependency: LaTeXStrings...")
        Pkg.add("LaTeXStrings")
    end
    println("✓ Environment is up to date with all packages installed.")
'

# 3. Directory Structure Initialization
echo "📁 Verification of local directory structures..."
mkdir -p "$DATA_DIR"
mkdir -p "$REPORTS_DIR"
mkdir -p "$DRAFTS_DIR"
echo "✓ Core output and asset directories confirmed."

# --- Operational Execution Routines ---
print_usage() {
    echo "Usage: $0 [flag]"
    echo "Flags:"
    echo "  --init          Install packages, initialize directories, and exit."
    echo "  --day YYMMDD    Run diagnostic pipeline for a specific campaign day (e.g., --day 991024)."
    echo "  --month         Execute batch mode sweep across all active October campaign records."
    echo "  --help          Display this help instruction matrix."
}

# Check for empty arguments
if [ -z "$1" ]; then
    print_usage
    exit 0
fi

case "$1" in
    --init)
        echo "✅ Initialization complete. Environment is ready for pipeline calculations."
        exit 0
        ;;

    --day)
        if [ -z "$2" ]; then
            echo "❌ Error: --day flag requires a target campaign date (e.g., $0 --day 991024)."
            exit 1
        fi
        TARGET_DAY="$2"
        echo "🚀 Launching spectral manifold diagnostics for target day: 19$TARGET_DAY..."

        # Executes the diagnostic pipeline script, passing the target day parameter directly to ARGS
        julia --project="$PROJECT_DIR" "$SCRIPT_DIR/run_diagnostics.jl" "$TARGET_DAY"
        echo "✅ Diagnostic execution successfully finished for day 19$TARGET_DAY."
        ;;

    --month)
        echo "🚀 Initiating full monthly verification run (October 1–31, 1999)..."
        # If your primary monthly loop batch execution is handled by another script:
        if [ -f "$SCRIPT_DIR/RunCampaignPipeline.jl" ]; then
            julia --project="$PROJECT_DIR" "$SCRIPT_DIR/RunCampaignPipeline.jl"
        else
            echo "⚠️ Warning: $SCRIPT_DIR/RunCampaignPipeline.jl not found."
            echo "Falling back to running your diagnostics utility as a standalone process..."
            julia --project="$PROJECT_DIR" "$SCRIPT_DIR/run_diagnostics.jl"
        fi
        echo "✅ Full month processing pipeline finalized."
        ;;

    --help)
        print_usage
        exit 0
        ;;

    *)
        echo "❌ Error: Unknown execution flag '$1'"
        print_usage
        exit 1
        ;;
esac

echo "=========================================================="
echo "🎉 Operation pipeline workflow completed successfully."
echo "=========================================================="
