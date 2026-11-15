"""Compatibility between upstream measurement workflows and fork transport."""

import numpy as np
import pandas as pd
import pytest
import iadpython as iad


def table(gs):
    mu = np.linspace(-1, 1, 401)
    return pd.DataFrame({i: (1-g*g)/(2*(1+g*g-2*g*mu)**1.5)
                         for i, g in enumerate(gs)}, index=mu)


@pytest.mark.parametrize("kind", ["TABULATED", "MOMENTS"])
def test_point_at_slices_custom_phase(kind):
    data = table([0.2, 0.6]) if kind == "TABULATED" else np.array([0.2, 0.6])**np.arange(33)[:, None]
    sample = iad.Sample(a=[0.8, 0.9], b=[0.4, 0.8], pf_type=kind, pf_data=data)
    experiment = iad.Experiment(sample=sample)
    expected = np.asarray(sample.rt())
    for i in range(2):
        point = experiment.point_at(i)
        np.testing.assert_allclose(np.asarray(point.sample.rt()).reshape(4), expected[:, i], rtol=1e-12)
    assert sample.pf_data.shape[1] == 2


def test_point_at_preserves_mismatched_spectral_layers():
    sample = iad.Sample(layers=[
        iad.Layer(a=[0.8, 0.9], b=0.4, n=[1.3, 1.4], pf_type="TABULATED", pf_data=table([0.2, 0.6])),
        iad.Layer(a=0.6, b=[0.2, 0.5], n=1.6, pf_type="MOMENTS", pf_data=0.3**np.arange(129)),
    ])
    expected = np.asarray(sample.rt())
    ballistic = np.asarray(sample.unscattered_rt())
    for i in range(2):
        point = iad.Experiment(sample=sample).point_at(i)
        np.testing.assert_allclose(np.asarray(point.sample.rt()).reshape(4), expected[:, i], rtol=1e-12)
        np.testing.assert_allclose(point.sample.unscattered_rt(), ballistic[:, i], rtol=1e-12)
    assert sample.layers[0].pf_data.shape[1] == 2


def test_measurement_uses_external_angle_and_absorbing_slides():
    sample = iad.Sample(a=0.8, b=0.7, n=1.4, n_above=1.5, n_below=1.6,
                        n_outer_above=1.2, n_outer_below=1.1)
    sample.nu_0 = 0.65
    sample.b_above, sample.b_below = 0.1, 0.2
    ru, tu = iad.specular_rt(1.5, 1.4, 1.6, 0.7, 0.65, b_top=0.1, b_bot=0.2,
                            n_outer_top=1.2, n_outer_bot=1.1, n_slab_transport=1.4)
    np.testing.assert_allclose(sample.unscattered_rt(), [ru, tu])
    experiment = iad.Experiment(sample=sample)
    experiment.fraction_of_rc_in_mr = experiment.fraction_of_tc_in_mt = 0
    np.testing.assert_allclose(experiment.measured_rt_from_raw(0.5, 0.6, 0.4, 0.3),
                               [0.5-ru, 0.6-tu])


@pytest.mark.parametrize("sample", [iad.Sample(a=0.8), iad.Sample(a=[0.7, 0.8]),
    iad.Sample(layers=[iad.Layer(a=0.8, n=1.3), iad.Layer(a=0.7, n=1.5)])])
def test_forward_counter_does_not_change_results(sample):
    first = sample.rt()
    assert sample.rt_evals == 1
    np.testing.assert_allclose(sample.rt(), first, rtol=0, atol=0)
    assert sample.rt_evals == 2


def test_legacy_hg_rejects_custom_phase():
    with pytest.raises(ValueError, match="phase_legendre"):
        iad.hg_legendre(iad.Sample(pf_type="TABULATED", pf_data=table([0.2])))


def test_matched_grazing_interface():
    assert iad.fresnel.interface_rt(1.4, 0, 1.4) == (0, 1)
    r, t = iad.fresnel.interface_rt(1.4, np.array([0, 0.5, 1]), 1.4)
    np.testing.assert_array_equal(r, 0)
    np.testing.assert_array_equal(t, 1)


def test_boundary_index_tracks_loader_updates_but_preserves_override():
    sample = iad.Sample()
    sample.n = 1.4
    reference = iad.Sample(n=1.4)
    np.testing.assert_allclose(sample.rt(), reference.rt())
    np.testing.assert_allclose(sample.unscattered_rt(), reference.unscattered_rt())
    sample.n_sample_boundary = 1.5 + 0.01j
    sample.n = 1.6
    assert sample.n_sample_boundary == 1.5 + 0.01j
    sample.n_sample_boundary = None
    assert sample.n_sample_boundary == 1.6


def test_mixed_phase_stacks_match_premerge_fork():
    import json
    from pathlib import Path
    baseline = json.loads((Path(__file__).parent / "data" / "fork_forward_reference.json").read_text())
    for outer, expected in zip([(1., 1.), (1.2, 1.1)], baseline["results"]):
        sample = iad.Sample(layers=[
            iad.Layer(a=[.8, .9], b=.4, n=[1.3, 1.4], pf_type="TABULATED", pf_data=table([.2, .6])),
            iad.Layer(a=.7, b=[.3, .5], n=1.6, pf_type="MOMENTS", pf_data=.3**np.arange(129)),
            iad.Layer(a=.4, b=.2, n=1.2, g=.1),
        ], n_outer_above=outer[0], n_outer_below=outer[1], quad_pts=8)
        np.testing.assert_allclose(sample.rt(), expected[0], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(sample.unscattered_rt(), expected[1], rtol=1e-12, atol=1e-14)
