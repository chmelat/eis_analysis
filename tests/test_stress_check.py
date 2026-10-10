"""The stress test's --check: case by case, but counted for the groups whose
failing cases depend on the platform (tests/stress.py, COUNTED)."""

from tests.stress import check_baseline


def _results(failures, noisy=()):
    """Fake run of cases rc/0..rc/39, noise-free except the indices in noisy,
    with the given 'rc/i analysis:invariant' failures."""
    results = [{'family': 'rc', 'index': i, 'noise': f'{0.01 if i in noisy else 0:g} proportional',
                'checks': []} for i in range(40)]
    for key in failures:
        case, check = key.split()
        analysis, invariant = check.split(':')
        results[int(case.split('/')[1])]['checks'].append(
            {'analysis': analysis, 'invariant': invariant, 'status': 'fail'})
    return results


BASE = {'failures': [f'rc/{i} linkk:B+' for i in range(9)] + ['rc/30 rinf:Irange'],
        'rates': {}}


def test_swapped_cases_of_a_counted_group_pass():
    swapped = [f'rc/{i} linkk:B+' for i in range(10, 19)] + ['rc/30 rinf:Irange']
    assert check_baseline(_results(swapped), False, BASE)


def test_count_above_the_allowance_fails():
    # 9 in the baseline: allowed 9 + max(sqrt(9), 3) = 12
    assert check_baseline(_results([f'rc/{i} linkk:B+' for i in range(12)]), False, BASE)
    assert not check_baseline(_results([f'rc/{i} linkk:B+' for i in range(13)]), False, BASE)


def test_empty_counted_group_allows_the_slack():
    fits = [f'rc/{i} fit:Cmix' for i in range(31, 35)]
    assert check_baseline(_results(BASE['failures'] + fits[:3]), False, BASE)
    assert not check_baseline(_results(BASE['failures'] + fits), False, BASE)


def test_linkk_on_a_noisy_case_is_checked_by_case():
    # Lin-KK fails B/C only on noise-free spectra (known limit 1)
    assert not check_baseline(_results(BASE['failures'] + ['rc/35 linkk:B+'], noisy=(35,)),
                              False, BASE)
    # a fit is counted at any noise (limit 3 has noisy ill-posed fits)
    assert check_baseline(_results(BASE['failures'] + ['rc/35 fit:B+'], noisy=(35,)), False, BASE)


def test_new_failure_outside_counted_groups_fails():
    assert not check_baseline(_results(BASE['failures'] + ['rc/31 rinf:Irange']), False, BASE)
    assert not check_baseline(_results(BASE['failures'] + ['rc/31 linkk:D']), False, BASE)
