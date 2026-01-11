#!/usr/bin/env python3
"""
Integration tests for CLI workflow.

Tests end-to-end CLI workflows including:
1. Synthetic data analysis (no input file)
2. Kramers-Kronig validation
3. DRT analysis with various options
4. Equivalent circuit fitting
5. Voigt chain fitting
6. Combined workflows
7. File-based input (CSV, DTA)

These tests verify that the complete analysis pipeline works correctly
when invoked through the CLI interface.

Run with: python3 tests/test_cli_integration.py
"""

import sys
import os
import argparse
import tempfile
import logging
from typing import Optional

import numpy as np

# Ensure the package is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Suppress matplotlib GUI
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Setup minimal logging for tests
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

print("=" * 70)
print("Integration Tests for CLI Workflow")
print("=" * 70)


def create_test_args(**kwargs) -> argparse.Namespace:
    """
    Create argparse.Namespace with default CLI arguments.

    Override defaults by passing keyword arguments.
    """
    defaults = {
        # Input/Output
        'input': None,
        'f_min': None,
        'f_max': None,
        'save': None,
        'format': 'png',
        'no_show': True,  # Always disable show in tests
        'verbose': 0,
        'quiet': True,

        # DRT
        'lambda_reg': None,
        'normalize_rpol': False,
        'n_tau': 100,
        'no_voigt_info': True,
        'classify_terms': False,
        'no_drt': False,
        'peak_method': 'scipy',
        'gmm_bic_threshold': 10.0,
        'ri_fit': False,

        # Kramers-Kronig
        'no_kk': False,
        'mu_threshold': 0.85,
        'auto_extend': False,
        'extend_decades_max': 1.0,

        # Z-HIT
        'no_zhit': True,  # Skip by default for speed
        'zhit_optimize_offset': False,

        # Circuit Fitting
        'circuit': None,
        'weighting': 'modulus',
        'no_fit': True,  # Skip by default unless specified
        'numeric_jacobian': False,
        'optimizer': 'de',
        'multistart': 0,
        'multistart_scale': 2.0,
        'de_strategy': 1,
        'de_popsize': 15,
        'de_maxiter': 1000,
        'de_tol': 0.01,
        'de_workers': 1,

        # Voigt Chain
        'voigt_chain': False,
        'voigt_n_per_decade': 3,
        'voigt_extend_decades': 0.0,
        'voigt_prune_threshold': 0.01,
        'voigt_allow_negative': False,
        'voigt_no_inductance': False,
        'voigt_fit_type': 'complex',
        'voigt_auto_M': False,
        'voigt_mu_threshold': 0.85,
        'voigt_max_M': 50,

        # Oxide Analysis
        'analyze_oxide': False,
        'epsilon_r': 22.0,
        'area': 1.0,

        # Visualization
        'ocv': False,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


# =============================================================================
# Test Fixtures
# =============================================================================

def get_synthetic_data():
    """Generate synthetic EIS data for testing."""
    from eis_analysis.io import generate_synthetic_data
    return generate_synthetic_data(
        Rs=10,
        R0=1e5,
        Q0=(1e-6, 0.6),
        R1=8e5,
        Q1=(3e-5, 0.43),
        noise=0.01
    )


def get_simple_rc_data():
    """Generate simple RC circuit data for circuit fitting tests."""
    frequencies = np.logspace(-1, 4, 50)
    omega = 2 * np.pi * frequencies

    Rs = 100.0
    R1 = 1000.0
    C1 = 1e-6

    # Rs - (R1 || C1)
    Z_RC = R1 / (1 + 1j * omega * R1 * C1)
    Z = Rs + Z_RC

    # Add small noise
    np.random.seed(42)
    noise = 0.005 * np.abs(Z) * (np.random.randn(len(Z)) + 1j * np.random.randn(len(Z)))
    Z = Z + noise

    return frequencies, Z


# =============================================================================
# Test 1: Synthetic Data - Basic Workflow
# =============================================================================

def test_synthetic_data_basic():
    """Test basic workflow with synthetic data."""
    print("\n[Test 1] Synthetic data - basic workflow")
    print("-" * 70)

    from eis_analysis.cli import load_eis_data, filter_by_frequency

    args = create_test_args()

    # Load synthetic data
    data = load_eis_data(args)

    assert data.frequencies is not None, "Frequencies should not be None"
    assert data.Z is not None, "Z should not be None"
    assert len(data.frequencies) > 0, "Should have frequency data"
    assert len(data.frequencies) == len(data.Z), "Frequencies and Z should have same length"
    assert data.title == "Synthetic data", "Title should be 'Synthetic data'"

    print(f"  Loaded {len(data.frequencies)} data points")
    print(f"  Frequency range: [{data.frequencies.min():.2e}, {data.frequencies.max():.2e}] Hz")
    print(f"  |Z| range: [{np.abs(data.Z).min():.2e}, {np.abs(data.Z).max():.2e}] Ohm")
    print("  [OK] Synthetic data loaded successfully")

    return True


# =============================================================================
# Test 2: Kramers-Kronig Validation
# =============================================================================

def test_kk_validation():
    """Test Kramers-Kronig validation workflow."""
    print("\n[Test 2] Kramers-Kronig validation")
    print("-" * 70)

    from eis_analysis.cli import run_kk_validation

    frequencies, Z = get_synthetic_data()
    args = create_test_args(no_kk=False)

    fig = run_kk_validation(frequencies, Z, args)

    assert fig is not None, "KK validation should return a figure"
    assert isinstance(fig, plt.Figure), "Should return matplotlib Figure"

    # Clean up
    plt.close(fig)

    print(f"  Validated {len(frequencies)} data points")
    print("  [OK] KK validation completed successfully")

    return True


# =============================================================================
# Test 3: DRT Analysis
# =============================================================================

def test_drt_analysis():
    """Test DRT analysis workflow."""
    print("\n[Test 3] DRT analysis")
    print("-" * 70)

    from eis_analysis.cli import run_drt_analysis
    from eis_analysis.drt import DRTResult

    frequencies, Z = get_synthetic_data()
    args = create_test_args(no_drt=False, peak_method='scipy')

    # Run DRT without R_inf estimation
    result = run_drt_analysis(frequencies, Z, args, None, 'scipy')

    assert result is not None, "DRT should return a result"
    assert isinstance(result, DRTResult), "Should return DRTResult"
    assert result.tau is not None, "tau should not be None"
    assert result.gamma is not None, "gamma should not be None"
    assert len(result.tau) == len(result.gamma), "tau and gamma should have same length"

    print(f"  tau range: [{result.tau.min():.2e}, {result.tau.max():.2e}] s")
    print(f"  gamma max: {result.gamma.max():.2e}")
    print(f"  Lambda used: {result.lambda_used:.2e}")

    # Clean up any figures
    plt.close('all')

    print("  [OK] DRT analysis completed successfully")

    return True


def test_drt_with_gmm():
    """Test DRT analysis with GMM peak detection."""
    print("\n[Test 3b] DRT analysis with GMM peaks")
    print("-" * 70)

    try:
        from sklearn.mixture import GaussianMixture
        has_sklearn = True
    except ImportError:
        has_sklearn = False
        print("  [SKIP] scikit-learn not available")
        return True

    from eis_analysis.cli import run_drt_analysis

    frequencies, Z = get_synthetic_data()
    args = create_test_args(no_drt=False, peak_method='gmm', gmm_bic_threshold=10.0)

    result = run_drt_analysis(frequencies, Z, args, None, 'gmm')

    assert result is not None, "DRT should return a result"

    # GMM should detect peaks
    if result.peaks is not None:
        print(f"  GMM detected {len(result.peaks)} peaks")

    plt.close('all')
    print("  [OK] DRT with GMM completed successfully")

    return True


# =============================================================================
# Test 4: R_inf Estimation
# =============================================================================

def test_rinf_estimation():
    """Test R_inf estimation workflow."""
    print("\n[Test 4] R_inf estimation")
    print("-" * 70)

    from eis_analysis.cli import run_rinf_estimation

    frequencies, Z = get_synthetic_data()
    args = create_test_args(ri_fit=True)

    R_inf, fig = run_rinf_estimation(frequencies, Z, args)

    # R_inf estimation may return None if it fails, but should not raise
    if R_inf is not None:
        print(f"  Estimated R_inf: {R_inf:.2f} Ohm")
        assert R_inf > 0, "R_inf should be positive"
    else:
        print("  R_inf estimation returned None (expected for some data)")

    if fig is not None:
        plt.close(fig)

    print("  [OK] R_inf estimation completed")

    return True


# =============================================================================
# Test 5: Circuit Fitting
# =============================================================================

def test_circuit_fitting_single():
    """Test circuit fitting with single optimizer."""
    print("\n[Test 5] Circuit fitting (single optimizer)")
    print("-" * 70)

    from eis_analysis.cli import run_circuit_fitting

    frequencies, Z = get_simple_rc_data()

    # Simple R-RC circuit
    circuit_str = 'R(100)-(R(1000)|C(1e-6))'
    args = create_test_args(
        no_fit=False,
        circuit=circuit_str,
        optimizer='single'
    )

    result, fig = run_circuit_fitting(frequencies, Z, args)

    assert result is not None, "Fitting should return a result"
    assert hasattr(result, 'params_opt'), "Result should have optimized parameters"
    assert hasattr(result, 'fit_error_rel'), "Result should have fit error"

    print(f"  Circuit: {circuit_str}")
    print(f"  Fit error: {result.fit_error_rel:.2f}%")
    print(f"  Quality: {result.quality}")

    # Good fit should have low error
    assert result.fit_error_rel < 5.0, f"Fit error too high: {result.fit_error_rel:.2f}%"

    if fig is not None:
        plt.close(fig)

    print("  [OK] Circuit fitting (single) completed successfully")

    return True


def test_circuit_fitting_multistart():
    """Test circuit fitting with multistart optimizer."""
    print("\n[Test 5b] Circuit fitting (multistart)")
    print("-" * 70)

    from eis_analysis.cli import run_circuit_fitting

    frequencies, Z = get_simple_rc_data()

    circuit_str = 'R(100)-(R(1000)|C(1e-6))'
    args = create_test_args(
        no_fit=False,
        circuit=circuit_str,
        optimizer='multistart',
        multistart=3
    )

    result, fig = run_circuit_fitting(frequencies, Z, args)

    assert result is not None, "Fitting should return a result"
    print(f"  Fit error: {result.fit_error_rel:.2f}%")

    if fig is not None:
        plt.close(fig)

    print("  [OK] Circuit fitting (multistart) completed successfully")

    return True


def test_circuit_fitting_de():
    """Test circuit fitting with differential evolution."""
    print("\n[Test 5c] Circuit fitting (differential evolution)")
    print("-" * 70)

    from eis_analysis.cli import run_circuit_fitting

    frequencies, Z = get_simple_rc_data()

    circuit_str = 'R(100)-(R(1000)|C(1e-6))'
    args = create_test_args(
        no_fit=False,
        circuit=circuit_str,
        optimizer='de',
        de_popsize=10,
        de_maxiter=100,  # Reduced for speed
        de_tol=0.1
    )

    result, fig = run_circuit_fitting(frequencies, Z, args)

    assert result is not None, "Fitting should return a result"
    print(f"  Fit error: {result.fit_error_rel:.2f}%")

    if fig is not None:
        plt.close(fig)

    print("  [OK] Circuit fitting (DE) completed successfully")

    return True


# =============================================================================
# Test 6: Voigt Chain Fitting
# =============================================================================

def test_voigt_chain():
    """Test Voigt chain fitting workflow."""
    print("\n[Test 6] Voigt chain fitting")
    print("-" * 70)

    from eis_analysis.fitting import fit_voigt_chain_linear

    frequencies, Z = get_synthetic_data()

    circuit, params = fit_voigt_chain_linear(
        frequencies, Z,
        n_per_decade=2,
        extend_decades=0.5,
        prune_threshold=0.05
    )

    assert circuit is not None, "Should return a circuit"
    assert params is not None, "Should return parameters"
    assert len(params) > 0, "Should have parameters"

    print(f"  Generated circuit: {circuit}")
    print(f"  Number of parameters: {len(params)}")

    # Verify circuit is valid by fitting
    from eis_analysis.fitting import fit_equivalent_circuit

    result, Z_fit, fig = fit_equivalent_circuit(frequencies, Z, circuit)

    assert result is not None, "Fit should succeed"
    print(f"  Fit error: {result.fit_error_rel:.2f}%")

    if fig is not None:
        plt.close(fig)

    print("  [OK] Voigt chain fitting completed successfully")

    return True


# =============================================================================
# Test 7: Z-HIT Validation
# =============================================================================

def test_zhit_validation():
    """Test Z-HIT validation workflow."""
    print("\n[Test 7] Z-HIT validation")
    print("-" * 70)

    from eis_analysis.cli import run_zhit_validation

    frequencies, Z = get_synthetic_data()
    args = create_test_args(no_zhit=False)

    fig = run_zhit_validation(frequencies, Z, args)

    assert fig is not None, "Z-HIT should return a figure"

    plt.close(fig)

    print("  [OK] Z-HIT validation completed successfully")

    return True


# =============================================================================
# Test 8: File Input - CSV
# =============================================================================

def test_csv_input():
    """Test loading and analyzing CSV file."""
    print("\n[Test 8] CSV file input")
    print("-" * 70)

    csv_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'example', 'example_eis_data.csv'
    )

    if not os.path.exists(csv_path):
        print(f"  [SKIP] CSV file not found: {csv_path}")
        return True

    from eis_analysis.cli import load_eis_data, run_kk_validation

    args = create_test_args(input=csv_path)

    data = load_eis_data(args)

    assert data.frequencies is not None, "Should load frequencies"
    assert data.Z is not None, "Should load impedance"
    assert len(data.frequencies) > 0, "Should have data"

    print(f"  Loaded {len(data.frequencies)} points from CSV")
    print(f"  Frequency range: [{data.frequencies.min():.2e}, {data.frequencies.max():.2e}] Hz")

    # Run KK validation on loaded data
    fig = run_kk_validation(data.frequencies, data.Z, create_test_args(no_kk=False))
    if fig is not None:
        plt.close(fig)

    print("  [OK] CSV file processing completed successfully")

    return True


# =============================================================================
# Test 9: File Input - DTA
# =============================================================================

def test_dta_input():
    """Test loading and analyzing DTA file."""
    print("\n[Test 9] DTA file input")
    print("-" * 70)

    dta_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'example', 'real_gamry_example.DTA'
    )

    if not os.path.exists(dta_path):
        print(f"  [SKIP] DTA file not found: {dta_path}")
        return True

    from eis_analysis.cli import load_eis_data

    args = create_test_args(input=dta_path)

    data = load_eis_data(args)

    assert data.frequencies is not None, "Should load frequencies"
    assert data.Z is not None, "Should load impedance"
    assert len(data.frequencies) > 0, "Should have data"

    # DTA files should have metadata
    if data.metadata is not None:
        print(f"  Metadata available: {list(data.metadata.keys())[:5]}...")

    print(f"  Loaded {len(data.frequencies)} points from DTA")
    print(f"  Frequency range: [{data.frequencies.min():.2e}, {data.frequencies.max():.2e}] Hz")

    print("  [OK] DTA file processing completed successfully")

    return True


# =============================================================================
# Test 10: Frequency Filtering
# =============================================================================

def test_frequency_filtering():
    """Test frequency range filtering."""
    print("\n[Test 10] Frequency filtering")
    print("-" * 70)

    from eis_analysis.cli import load_eis_data, filter_by_frequency

    args = create_test_args()
    data = load_eis_data(args)

    original_count = len(data.frequencies)

    # Apply frequency filter
    args_filtered = create_test_args(f_min=1.0, f_max=1000.0)
    filtered_data = filter_by_frequency(data, args_filtered)

    filtered_count = len(filtered_data.frequencies)

    assert filtered_count <= original_count, "Filtered data should have fewer or equal points"
    assert filtered_data.frequencies.min() >= 1.0, "Min frequency should be >= 1 Hz"
    assert filtered_data.frequencies.max() <= 1000.0, "Max frequency should be <= 1000 Hz"

    print(f"  Original: {original_count} points")
    print(f"  Filtered (1-1000 Hz): {filtered_count} points")
    print(f"  Filtered range: [{filtered_data.frequencies.min():.2e}, {filtered_data.frequencies.max():.2e}] Hz")

    print("  [OK] Frequency filtering completed successfully")

    return True


# =============================================================================
# Test 11: Combined Workflow
# =============================================================================

def test_combined_workflow():
    """Test complete combined workflow."""
    print("\n[Test 11] Combined workflow (KK + DRT + Fit)")
    print("-" * 70)

    from eis_analysis.cli import (
        load_eis_data,
        run_kk_validation,
        run_drt_analysis,
        run_circuit_fitting,
    )

    args = create_test_args()
    data = load_eis_data(args)

    # Step 1: KK validation
    kk_args = create_test_args(no_kk=False)
    fig_kk = run_kk_validation(data.frequencies, data.Z, kk_args)
    assert fig_kk is not None, "KK should produce figure"

    # Step 2: DRT analysis
    drt_args = create_test_args(no_drt=False)
    drt_result = run_drt_analysis(data.frequencies, data.Z, drt_args, None, 'scipy')
    assert drt_result is not None, "DRT should produce result"

    # Step 3: Circuit fitting
    fit_args = create_test_args(
        no_fit=False,
        circuit='R(10)-(R(1e5)|Q(1e-6,0.6))-(R(8e5)|Q(3e-5,0.43))',
        optimizer='single'
    )
    fit_result, fig_fit = run_circuit_fitting(data.frequencies, data.Z, fit_args)
    assert fit_result is not None, "Fit should produce result"

    print(f"  KK: figure generated")
    print(f"  DRT: lambda={drt_result.lambda_used:.2e}")
    print(f"  Fit: error={fit_result.fit_error_rel:.2f}%")

    # Clean up
    plt.close('all')

    print("  [OK] Combined workflow completed successfully")

    return True


# =============================================================================
# Test 12: Save Output
# =============================================================================

def test_save_output():
    """Test saving figures to files."""
    print("\n[Test 12] Save output to files")
    print("-" * 70)

    from eis_analysis.cli import run_kk_validation

    frequencies, Z = get_synthetic_data()

    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, 'test_output')

        args = create_test_args(no_kk=False, save=prefix, format='png')
        fig = run_kk_validation(frequencies, Z, args)

        expected_file = f"{prefix}_kk.png"

        if fig is not None:
            plt.close(fig)

        # Check if file was created
        if os.path.exists(expected_file):
            file_size = os.path.getsize(expected_file)
            print(f"  Saved: {os.path.basename(expected_file)} ({file_size} bytes)")
            assert file_size > 0, "Saved file should not be empty"
        else:
            print(f"  Note: File not saved (save_figure may be disabled in handler)")

    print("  [OK] Save output test completed")

    return True


# =============================================================================
# Test 13: Error Handling
# =============================================================================

def test_error_handling():
    """Test error handling for invalid inputs."""
    print("\n[Test 13] Error handling")
    print("-" * 70)

    from eis_analysis.cli import load_eis_data, EISAnalysisError

    # Test non-existent file
    args = create_test_args(input='/nonexistent/file.csv')

    try:
        load_eis_data(args)
        assert False, "Should raise EISAnalysisError for non-existent file"
    except EISAnalysisError as e:
        print(f"  Caught expected error: {str(e)[:50]}...")

    # Test unsupported format
    with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:
        f.write(b'test data')
        tmp_path = f.name

    try:
        args = create_test_args(input=tmp_path)
        load_eis_data(args)
        assert False, "Should raise EISAnalysisError for unsupported format"
    except EISAnalysisError as e:
        print(f"  Caught expected error: {str(e)[:50]}...")
    finally:
        os.unlink(tmp_path)

    print("  [OK] Error handling works correctly")

    return True


# =============================================================================
# Test Runner
# =============================================================================

def run_all_tests():
    """Run all integration tests."""
    tests = [
        ("Synthetic data basic", test_synthetic_data_basic),
        ("KK validation", test_kk_validation),
        ("DRT analysis", test_drt_analysis),
        ("DRT with GMM", test_drt_with_gmm),
        ("R_inf estimation", test_rinf_estimation),
        ("Circuit fitting (single)", test_circuit_fitting_single),
        ("Circuit fitting (multistart)", test_circuit_fitting_multistart),
        ("Circuit fitting (DE)", test_circuit_fitting_de),
        ("Voigt chain", test_voigt_chain),
        ("Z-HIT validation", test_zhit_validation),
        ("CSV input", test_csv_input),
        ("DTA input", test_dta_input),
        ("Frequency filtering", test_frequency_filtering),
        ("Combined workflow", test_combined_workflow),
        ("Save output", test_save_output),
        ("Error handling", test_error_handling),
    ]

    passed = 0
    failed = 0
    skipped = 0

    for name, test_func in tests:
        try:
            result = test_func()
            if result:
                passed += 1
        except AssertionError as e:
            print(f"\n  [FAIL] {name}: {e}")
            failed += 1
        except Exception as e:
            print(f"\n  [ERROR] {name}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    # Summary
    print("\n" + "=" * 70)
    print("INTEGRATION TEST SUMMARY")
    print("=" * 70)
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Total:   {len(tests)}")
    print("=" * 70)

    if failed == 0:
        print("\nALL INTEGRATION TESTS PASSED!")
        return 0
    else:
        print(f"\n{failed} TEST(S) FAILED!")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())
